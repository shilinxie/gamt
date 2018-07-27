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


#' Best model
#'
#' @param x numeric vector of predictor variable
#' @param y numeric vector of response variable
#' @param ... other arguments to pass to the mgcv::gam()
#'
#' @return
#' @export
#'
#' @examples
#' x <- seq(0, 50, 1)
#' y <- ((runif(1, 10, 20) * x) / (runif(1, 0, 10) + x)) + rnorm(51, 0, 1)
#' ## fit the nonlinear model using maximum likelihood ("ML")
#' bb <- best_model(x = x, y = y, method = "ML")
#' summary(bb)

best_model <- function(x, y, ...){

  ## Check data
  if(length(y) != length(x)){
    stop("x and y must be of equal length")
  }

  ## add default model arguments if not otherwise supplied
  if(!exists("bs")) {
    bs <- "tp"
  }

  if(!exists("k")) {
    k <- 4
  }

  if(!exists("method")) {
    method = "GCV.Cp"
  }

  try_model <- function(model_type) {
    out <- tryCatch(
      {
        ## Define random component for ac models
        rand <- rep(1, length(x))

        ## Define model formula
        gam_formula <- as.formula(sprintf("%s ~ s(%s, bs = '%s', k = %s)",
                                          "y",
                                          "x",
                                          bs,
                                          k))

        lm_formula <- as.formula(sprintf("%s ~ %s",
                                         "y",
                                         "x"))
        ## run the models
        switch(model_type,
               gam = {mgcv::gam(gam_formula,
                                method = method,
                                se = T)
               },
               lm = {mgcv::gam(lm_formula,
                               method = method,
                               se = T)
               },
               gamm = {mgcv::gamm(gam_formula,
                                  se = T,
                                  correlation = nlme::corAR1(form = ~ 1))
               },
               lmac = {mgcv::gamm(lm_formula,
                                  random = list(rand = ~ 1),
                                  se = T,
                                  correlation = nlme::corAR1(form = ~ 1))
               },
               stop("Enter a valid model_type!")
        )
      },
      error = function(cond) {
        message(sprintf("The %s caused an error!\n%s\n",
                        model_type,
                        cond))
        return(NA)
      },
      warning = function(cond) {
        message(sprintf("The %s caused a warning!\n%s\n",
                        model_type,
                        cond))
        return(NA)
      }
    )
    return(out)
  }

  ## Run models
  gam_model   <- try_model(model_type = "gam")
  lm_model    <- try_model(model_type = "lm")
  gamm_model  <- try_model(model_type = "gamm")
  lmac_model  <- try_model(model_type = "lmac")

  ## Collect output
  gam_summary <- data.frame(edf = summary(gam_model)$edf,
                            d_gcv = summary(gam_model)$sp.criterion - summary(lm_model)$sp.criterion,
                            d_aic = abs(AICc(gam_model) - AICc(lm_model)),
                            stringsAsFactors = FALSE)

  gamm_summary <- data.frame(edf = summary(gamm_model$gam)$edf,
                            d_aic =  abs(summary(gamm_model$lme)$AIC-summary(lmac_model$lme)$AIC),
                            stringsAsFactors = FALSE)

  best_model <- NULL
  if(ac_prob(gam_model) <= 0.05) {  # probability of autocorrelation <= 0.05
    if(gamm_summary$edf >= 2.0 &    # EDF gamm > 2.0
       gamm_summary$d_aic >= 2.0) { # dAIC > 2.0
      best_model <- get("gamm_model")
      best_model$best_model <- "gamm"
    }
    if(gamm_summary$d_aic < 2.0) {
      best_model <- get("lmac_model")
      best_model$best_model <- "lmac"
    }
  }
  if(ac_prob(gam_model) > 0.05){    # probability of autocorrelation > 0.05
    if(gam_summary$edf >= 2.0 &     # EDF gam > 2.0
       gam_summary$d_gcv < 0 &      # dGCV negative
       gam_summary$d_aic >= 2.0) {  # dAIC > 2.0
      best_model <- get("gam_model")
      best_model$best_model <- "gam"
    }
    if(gam_summary$d_aic < 2.0) {
      best_model <- get("lm_model")
      best_model$best_model <- "lm"
    }
  }

  if(is.null(best_model)){
    stop("Something went wrong and the best model couldn't be identified.")
  }
  if(!is.null(best_model)){
    return(best_model)
  }
}


gamt <- function(formula, bootstraps, ...){


}
