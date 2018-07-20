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

best_model <- function(data){

  ## Check data

  all_summary <- data.frame("MODEL" = c("gam_model", "lm_model", "gamm_model", "lmac_model"),
                            stringsAsFactors = FALSE)
  gam_summary <- data.frame() %>%
    mutate(model = "gam",
           aicc = AICc(gam_model),
           loglik = logLik(gam_model),
           dev_expl = summary(gam_model)$dev.expl,
           edf = summary(gam_model)$edf,
           p_val = summary(gam_model)$s.pv,
           rsq = summary(gam_model)$r.sq,
           gcv = summary(gam_model)$sp.criterion,
           d_gcv = delta.GCV.gam.lm,
           d_aic = delta.AIC.gam.lm,
           d_dev_expl = dev.diff.gam.lm)


  ## Run models
  gam_model  <- mgcv::gam(age ~ s(circumference,
                                  bs = "tp",
                                  k = 4),
                          data = Orange,
                          optimMmethod = "GCV.Cp",
                          se = T)
  c("AICc", "logLik","dev.expl", "edf", "Pvalue","R-squared","GCV","delta.GCV","p.ac","delta.AIC","diff.dev.expl")

  lm_model  <- mgcv::gam(age ~ circumference,
                         data = Orange,
                         optimMmethod = "GCV.Cp",
                         se = T)

  gamm_model <- mgcv::gamm(age ~ s(circumference,
                                   bs = "tp",
                                   k = 4),
                           data = Orange,
                           optimMmethod = "GCV.Cp",
                           se = T,
                           correlation = nlme::corAR1(form = ~ts_len))

  lmac_model <-mgcv:: gamm(age ~ circumference,
                           random = list(rand = ~1),
                           data = Orange,
                           se = T,
                           correlation = nlme::corAR1(form = ~ts_len))

  ## Collect output




  ## Select best model

}


gamt <- function(formula, bootstraps, ...){

  gam_model <- mgcv::gam(formula, ...)


  gam1  <- gam(myresponse ~ s(mydriver, bs= "tp",k = ks), optimMmethod="GCV.Cp",se = T)
  linear <- gam(myresponse ~ mydriver, method = "GCV.Cp", se = T)
  dev.resid <- residuals(gam1,type='deviance')
  lme1 <- lme(dev.resid~1,random=~1|rand,correlation=corAR1(form=~year),method='ML')
  lm1 <- lm(dev.resid~1)
  p.ac <- 1-pchisq(2*(logLik(lme1)[1]-logLik(lm1)[1]),2)
  delta.GCV.gam.lm <- summary(gam1)$sp.criterion - summary(linear)$sp.criterion   #A negative value means the GAM with a smoother is a better model than the linear model
  delta.AIC.gam.lm <- AICc(gam1) - AICc(linear)                                   #A negative value means the GAM with a smoother is a better model than the linear model
  dev.diff.gam.lm <- summary(gam1)$dev.expl-summary(linear)$dev.expl


}
