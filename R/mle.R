#' Negative Log Likelihood
#'
#' Calculates the negative log-likelihood of a model for the data
#'
#' @param param a vector of values for the model paramters in the following order c(mux, sigmax, sig, alpha, b, beta_C)
#' @param LOD_x a single limit of detection for the censored covariate or a vector of LODs (specific to each observation)
#' @param LOD_y a single limit of detection for the outcome or a vector of LODs (specific to each observation)
#' @param x a vector of potentially left censored predictor values
#' @param y a vecotr of potentially left censored outcome values
#' @param censor_x a vector the length of x; 0 if x is observed, 1 if censored
#' @param censor_y a vector the length of y; 0 if y is observed, 1 if censored
#' @param C optional matrix or data frame of other covariates (not subject to LODs) included in the model
#' @param single_LODs do observations have varying limits of detection.  Logical.
#' @param x.dist distribution for the X covariate that is subject to a limit of detection.  'normal', 'gaussian' or 'lognormal'.
#'
#' @return negative log-likelihood of a model for the data
#'
#' @examples
#' x <- rnorm(100, mean = 1)
#' y <- 3*x + rnorm(100, mean = 0, sd=0.2)
#' censor_x <- ifelse(x < 0.2, 1, 0)
#' censor_y <- ifelse(y < 0, 1, 0)
#' NNL(param = c(1, 1, 0.2, alpha = 0, b=3), LOD_x = rep(0.2,100), LOD_y = rep(1,100), x = x, y = y, censor_x = censor_x, censor_y = censor_y, C=NULL, single_LODs=T, x.dist = 'normal')
#'
#' @export

NNL <- function(param, LOD_x, LOD_y, x, y, censor_x, censor_y, C=NULL, single_LODs=F, x.dist){
  # Calculate the likelihood
  if(x.dist=='lognormal'){
    if(is.null(C)==FALSE){
      LL <- rep(NA, length(x))
      if(length(LL[censor_x==0 & censor_y==0])>0) LL[censor_x==0 & censor_y==0] <-  LL1_lognormal(X=x[censor_x==0 & censor_y==0], Y=y[censor_x==0 & censor_y==0], mux = param[1], sigmax = param[2], sig=param[3], alpha=param[4], b=param[5], C[censor_x==0 & censor_y==0,], beta_C=param[6:length(param)])

      if(length(LL[censor_x==1 & censor_y==0])>0) LL[censor_x==1 & censor_y==0] <-  LL2_lognormal(Y=y[censor_x==1 & censor_y==0], LOD_X = LOD_x[censor_x==1 & censor_y==0], mux = param[1], sigmax = param[2], sig=param[3], alpha=param[4], b=param[5], C=C[censor_x==1 & censor_y==0,], beta_C=param[6:length(param)])

      if(length(LL[censor_x==0 & censor_y==1])>0) LL[censor_x==0 & censor_y==1] <-  LL3_lognormal(X=x[censor_x==0 & censor_y==1], LOD_Y = LOD_y[censor_x==0 & censor_y==1], mux = param[1], sigmax = param[2], sig=param[3], alpha=param[4], b=param[5], C=C[censor_x==0 & censor_y==1,], beta_C=param[6:length(param)])

      if(length(LL[censor_x==1 & censor_y==1])>0) LL[censor_x==1 & censor_y==1] <-  LL4_lognormal(LOD_X = LOD_x[censor_x==1 & censor_y==1], LOD_Y = LOD_y[censor_x==1 & censor_y==1], mux = param[1], sigmax = param[2], sig=param[3], alpha=param[4], b=param[5], C=C[censor_x==1 & censor_y==1,], beta_C=param[6:length(param)], single_LODs=F)

    }else{
      LL <- rep(NA, length(x))
      if(length(LL[censor_x==0 & censor_y==0])>0) LL[censor_x==0 & censor_y==0] <-  LL1_lognormal(X=x[censor_x==0 & censor_y==0], Y=y[censor_x==0 & censor_y==0], mux = param[1], sigmax = param[2], sig=param[3], alpha=param[4], b=param[5], C=NULL, beta_C=NULL)

      if(length(LL[censor_x==1 & censor_y==0])>0) LL[censor_x==1 & censor_y==0] <-  LL2_lognormal(Y=y[censor_x==1 & censor_y==0], LOD_X = LOD_x[censor_x==1 & censor_y==0], mux = param[1], sigmax = param[2], sig=param[3], alpha=param[4], b=param[5], C=NULL, beta_C=NULL)

      if(length(LL[censor_x==0 & censor_y==1])>0) LL[censor_x==0 & censor_y==1] <-  LL3_lognormal(X=x[censor_x==0 & censor_y==1], LOD_Y = LOD_y[censor_x==0 & censor_y==1], mux = param[1], sigmax = param[2], sig=param[3], alpha=param[4], b=param[5], C=NULL, beta_C=NULL)

      if(length(LL[censor_x==1 & censor_y==1]) >0 ) LL[censor_x==1 & censor_y==1] <-  LL4_lognormal(LOD_X = LOD_x[censor_x==1 & censor_y==1], LOD_Y = LOD_y[censor_x==1 & censor_y==1], mux = param[1], sigmax = param[2], sig=param[3], alpha=param[4], b=param[5], C=NULL, beta_C=NULL, single_LODs=single_LODs)
    }
  }else{

    if(is.null(C)==FALSE){
      LL <- rep(NA, length(x))
      if(length(LL[censor_x==0 & censor_y==0])>0) LL[censor_x==0 & censor_y==0] <-  LL1_normal(X=x[censor_x==0 & censor_y==0], Y=y[censor_x==0 & censor_y==0], mux = param[1], sigmax = param[2], sig=param[3], alpha=param[4], b=param[5], C[censor_x==0 & censor_y==0,], beta_C=param[6:length(param)])
      if(length(LL[censor_x==1 & censor_y==0])>0) LL[censor_x==1 & censor_y==0] <-  LL2_normal(Y=y[censor_x==1 & censor_y==0], LOD_X = LOD_x[censor_x==1 & censor_y==0], mux = param[1], sigmax = param[2], sig=param[3], alpha=param[4], b=param[5], C=C[censor_x==1 & censor_y==0,], beta_C=param[6:length(param)])
      if(length(LL[censor_x==0 & censor_y==1])>0) LL[censor_x==0 & censor_y==1] <-  LL3_normal(X=x[censor_x==0 & censor_y==1], LOD_Y = LOD_y[censor_x==0 & censor_y==1], mux = param[1], sigmax = param[2], sig=param[3], alpha=param[4], b=param[5], C=C[censor_x==0 & censor_y==1,], beta_C=param[6:length(param)])
      if(length(LL[censor_x==1 & censor_y==1])>0) LL[censor_x==1 & censor_y==1] <-  LL4_normal(LOD_X = LOD_x[censor_x==1 & censor_y==1], LOD_Y = LOD_y[censor_x==1 & censor_y==1], mux = param[1], sigmax = param[2], sig=param[3], alpha=param[4], b=param[5], C=C[censor_x==1 & censor_y==1,], beta_C=param[6:length(param)], single_LODs=F)
    }else{
      LL <- rep(NA, length(x))
      if(length(LL[censor_x==0 & censor_y==0])>0) LL[censor_x==0 & censor_y==0] <-  LL1_normal(X=x[censor_x==0 & censor_y==0], Y=y[censor_x==0 & censor_y==0], mux = param[1], sigmax = param[2], sig=param[3], alpha=param[4], b=param[5], C=NULL, beta_C=NULL)
      if(length(LL[censor_x==1 & censor_y==0])>0) LL[censor_x==1 & censor_y==0] <-  LL2_normal(Y=y[censor_x==1 & censor_y==0], LOD_X = LOD_x[censor_x==1 & censor_y==0], mux = param[1], sigmax = param[2], sig=param[3], alpha=param[4], b=param[5], C=NULL, beta_C=NULL)
      if(length(LL[censor_x==0 & censor_y==1])>0) LL[censor_x==0 & censor_y==1] <-  LL3_normal(X=x[censor_x==0 & censor_y==1], LOD_Y = LOD_y[censor_x==0 & censor_y==1], mux = param[1], sigmax = param[2], sig=param[3], alpha=param[4], b=param[5], C=NULL, beta_C=NULL)
      if(length(LL[censor_x==1 & censor_y==1])>0) LL[censor_x==1 & censor_y==1] <-  LL4_normal(LOD_X = LOD_x[censor_x==1 & censor_y==1], LOD_Y = LOD_y[censor_x==1 & censor_y==1], mux = param[1], sigmax = param[2], sig=param[3], alpha=param[4], b=param[5], C=NULL, beta_C=NULL, single_LODs=single_LODs)
    }

  }
  #MLE function will want to take the NEGATIVE log likelihood
  -sum(LL)
}


# MLE Fitting Function
#'
#' Maximum likelihood estimation of model paramters
#'
#' @param minuslogl a negative log likelihood function to minimize
#' @param start a vector of starting values for the parameters
#' @param method optimization method to use. See optim. Default is "Nelder-Mead".
#' @param nobs optional integer: the number of observations, to be used for e.g. computing BIC.
#' @param maxit maximum number of iterations
#' @param reltol relative tolerance for stopping optimization.  See optim.
#' @param ...	 Further arguments to pass to optim.
#'
#' @importClassesFrom stats4 mle
#'
#' @return an object of class mle
#'
#' @examples
#' x <- rnorm(100, mean = 1)
#' y <- 3*x + rnorm(100, mean = 0, sd=0.2)
#' censor_x <- ifelse(x < 0.2, 1, 0)
#' censor_y <- ifelse(y < 0, 1, 0)
#' mle2(start = c(mean(x), sd(x), sd(y), 0, 0), minuslogl = function(param){NNL(param, LOD_x = rep(0.2,100), LOD_y = rep(1,100), x = x, y = y, censor_x = censor_x, censor_y = censor_y, C=NULL, single_LODs=T, x.dist = 'normal')}, maxit = 1000)
#'
#'
#' @export
mle2 <- function (minuslogl, start = NULL, method = 'Nelder-Mead',
                  nobs, maxit=1000, reltol=1e-8, ...) {
  call <- match.call()

  oout <- optim(start, minuslogl, method = method, hessian = TRUE, control = list(maxit = maxit, reltol = reltol))

  coef <- oout$par
  vcov <- if (length(coef))
    solve(oout$hessian)
  else matrix(numeric(), 0L, 0L)
  min <- oout$value

  new("mle", call = call, coef = coef, fullcoef = coef,
      vcov = vcov, min = min, details = oout, minuslogl = minuslogl,
      nobs = if (missing(nobs))
        NA_integer_
      else nobs, method = method)
}





