#' AICc
#'
#' @param mod
#'
#' @return
#' @export
#'
#' @examples
AICc <- function(mod) {
  kc <- mod$rank
  nc <- length(mod$residuals)
  AICc <- round(mod$aic + (2 * kc *(kc + 1)/(nc - kc - 1)),
                3)
  return(AICc)
}


#' Deviance Explained
#'
#' @param mod
#'
#' @return
#' @export
#'
#' @examples
deviance_explained <- function(mod) {

  null_model  <- MASS::glmmPQL(age ~ 1,
                               random = list(rand = ~1),
                               data = Orange,
                               verbose = FALSE,
                               family='gaussian')
  dr <- sum(residuals(mod$gam)^2)
  dn0 <- sum(residuals(null_model)^2)

  return((dn0 - dr)/dn0)
}


#' Autocorrelation probability
#'
#' @param mod
#'
#' @return
#' @export
#'
#' @examples

ac_prob <- function(mod){

  if(!"gam" %in% class(mod)){
    stop("mod should be a gamObject")
  }

  ts_len <- seq_len(length(mod$y))
  rand <- rep(1, length(mod$y))

  deviance_resid <- residuals(mod,
                              type = 'deviance')

  lmac_resid <- nlme::lme(deviance_resid ~ 1,
                          random = ~1|rand,
                          correlation = nlme::corAR1(form = ~ts_len),
                          method = 'ML')

  lm_resid <- lm(deviance_resid ~ 1)
  ac_prob <- 1 - pchisq(2 * (logLik(lmac_resid)[1] - logLik(lm_resid)[1]),
                        2)
  return(ac_prob)
}


data(Orange)
Orange <- data.frame(Orange)
Orange$time <- seq(1:nrow(Orange))
ts_len <- Orange$time
rand <- rep(1, nrow(Orange))
data = Orange
response_name = "age"
term_name = "circumference"

best_model <- function(response_name, term_name, data, ...){

  ## Check data
  if(!"data.frame" %in% class(data)){
    stop("data must be supplied as a data.frame")
  }
  if(!all(c(response_name,
           term_name) %in% colnames(data))) {
    stop("Check that respone_name and term_name are both columns in dataframe")
  }
  if(!exists("bs")) {
    bs <- "tp"
  }
  if(!exists("k")) {
    k <- 4
  }
  if(!exists("optimMmethod")) {
    optimMmethod = "GCV.Cp"
  }

  rand <- rep(1, nrow(data))

  ## Define models
  gam_formula <- as.formula(sprintf("%s ~ s(%s, bs = '%s', k = %s)",
                                    response_name,
                                    term_name,
                                    bs,
                                    k))

  lm_formula <- as.formula(sprintf("%s ~ %s",
                                   response_name,
                                   term_name))

  ## Run models
  gam_model  <- mgcv::gam(gam_formula,
                          data = data,
                          optimMmethod = optimMmethod,
                          se = T)

  lm_model  <- mgcv::gam(lm_formula,
                         data = data,
                         optimMmethod = optimMmethod,
                         se = T)

  gamm_model <- mgcv::gamm(gam_formula,
                           data = data,
                           optimMmethod = optimMmethod,
                           se = T,
                           correlation = nlme::corAR1(form = ~ 1))

  lmac_model <-mgcv::gamm(lm_formula,
                          random = list(rand = ~ 1),
                          data = data,
                          se = T,
                          correlation = nlme::corAR1(form = ~ 1))

  ## Collect output
  gam_summary <- data.frame(#model = "gam_model",
                            # aicc = AICc(gam_model),
                            # loglik = logLik(gam_model),
                            # dev_expl = summary(gam_model)$dev.expl,
                            edf = summary(gam_model)$edf,
                            # p_val = summary(gam_model)$s.pv,
                            # rsq = summary(gam_model)$r.sq,
                            # gcv = summary(gam_model)$sp.criterion,
                            d_gcv = summary(gam_model)$sp.criterion - summary(lm_model)$sp.criterion ,
                            d_aic = AICc(gam_model) - AICc(lm_model) ,
                            # d_dev_expl = summary(gam_model)$dev.expl-summary(lm_model)$dev.expl,
                            stringsAsFactors = FALSE)

  gamm_summary <- data.frame(#model = "gamm_model",
                            # aicc = summary(gamm_model$lme)$AIC,
                            # loglik = summary(gamm_model$lme)$logLik,
                            # dev_expl = deviance_explained(gamm_model),
                            edf = summary(gamm_model$gam)$edf,
                            # p_val = summary(gamm_model$gam)$s.pv,
                            # rsq = summary(gamm_model$gam)$r.sq,
                            # gcv = NA,
                            # d_gcv = NA,
                            d_aic =  summary(gamm_model$lme)$AIC-summary(lmac_model$lme)$AIC,
                            # d_dev_expl = deviance_explained(gamm_model)-deviance_explained(lmac_model),
                            stringsAsFactors = FALSE)

  best_model <- NULL
  if(ac_prob(gam_model) <= 0.05) {  # probability of autocorrelation <= 0.05
    if(gamm_summary$edf >= 2.0 &    # EDF gamm > 2.0
       gamm_summary$d_aic >= 2.0) { # dAIC > 2.0
      best_model <- get("gamm_model")
    }
    if(gamm_summary$d_aic < 2.0) {
      best_model <- get("lmac_model")
    }
  }
  if(ac_prob(gam_model) > 0.05){    # probability of autocorrelation > 0.05
    if(gam_summary$edf >= 2.0 &     # EDF gam > 2.0
       gam_summary$d_gcv < 0 &      # dGCV negative
       gam_summary$d_aic >= 2.0) {  # dAIC > 2.0
      best_model <- get("gam_model")
    }
    if(gam_summary$d_aic < 2.0) {
      best_model <- get("lm_model")
    }
  }

  if(is.null(best_model)){
    stop("Something went wrong and the best model couldn't be identified.")
  }
  if(!is.null(best_model)){
    return(best_model)
  }
}


# lm_summary <- data.frame(model = "lm_model",
#                           aicc = AICc(lm_model),
#                           loglik = logLik(lm_model),
#                           dev_expl = summary(lm_model)$dev.expl,
#                           edf = summary(lm_model)$residual.df,
#                           p_val = summary(lm_model)$p.pv[[2]],
#                           rsq = summary(lm_model)$r.sq,
#                           gcv = summary(lm_model)$sp.criterion,
#                           d_gcv = NA,
#                           d_aic = NA,
#                           d_dev_expl = NA,
#                           stringsAsFactors = FALSE)

# lmac_summary <- data.frame(model = "lmac_model",
#                            aicc = summary(lmac_model$lme)$AIC,
#                            loglik = summary(lmac_model$lme)$logLik,
#                            dev_expl = deviance_explained(lmac_model),
#                            edf = summary(lmac_model$gam)$residual.df,
#                            p_val = summary(lmac_model$gam)$p.pv[[2]],
#                            rsq = summary(lmac_model$gam)$r.sq,
#                            gcv = NA,
#                            d_gcv = NA,
#                            d_aic = NA,
#                            d_dev_expl = NA,
#                            stringsAsFactors = FALSE)



gamt <- function(formula, bootstraps, ...){


}
