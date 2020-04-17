#' Log likelihood for a normally distributed covariate without covariate or outcome censoring
#'
#' Calculates the log-likelihood contribution for a normally distriubted covariate when neither the covariate or outcome are censored.
#'
#' @param X a normally distributed covariate
#' @param Y outcome
#' @param sig the standard deviation of the residuals
#' @param mux mean of the normally distributed covariate
#' @param sigmax standard deviation of the normally distributed covariate
#' @param alpha intercept: the average value of Y when all covariates = 0
#' @param b slope: the average change in Y for a 1 unit change in X
#' @param C optional vector of other covariates (not subject to LODs) included in the model
#' @param beta_C regression coefficients associated with C
#'
#' @return log-likelihood
#'
#' @examples
#' LL1_normal(X=2, Y=3, sig = 0.5, mux = 1, sigmax = 1, alpha = 0, b = 1)
#'
#' @export
LL1_normal <- function(X, Y, mux, sigmax, sig, alpha, b, C=NULL, beta_C=NULL){

  if(is.null(C)==F){e <- Y - alpha - b*X - c(as.matrix(C)%*%beta_C)
  }else{e <- Y - alpha - b*X}

  return((-(e^2)/(2*(sig^2))) - (1/2)*log((sig^2)*(sigmax^2)) - ((X - mux)^2)/((2*(sigmax^2))))

}


#' Log likelihood for a normally distributed left-censored covariate without outcome censoring
#'
#' Calculates the log-likelihood contribution for a normally distriubted covariate when the covariate is left censored and the outcome is not censored.
#'
#' @param LOD_X the limit of detection for the normally distributed covariate
#' @param Y outcome
#' @param sig the standard deviation of the residuals
#' @param mux mean of the normally distributed covariate
#' @param sigmax standard deviation of the normally distributed covariate
#' @param alpha intercept: the average value of Y when all covariates = 0
#' @param b slope: the average change in Y for a 1 unit change in X
#' @param C optional vector of other covariates (not subject to LODs) included in the model
#' @param beta_C regression coefficients associated with C
#'
#' @return log-likelihood
#'
#' @examples
#' LL2_normal(Y=3, LOD_X = 0.2, sig = 0.5, mux = 1, sigmax = 1, alpha = 0, b = 1)
#'
#' @export
LL2_normal <- function(Y, LOD_X, mux, sigmax, sig, alpha, b, C=NULL, beta_C=NULL){
  Q = sqrt((1/(sigmax^2)) + ((b^2)/(sig^2)))

  if(is.null(C)==F){emu <- Y - alpha - b*mux - c(as.matrix(C)%*%beta_C)
  }else{emu <- Y - alpha - b*mux}

  return((1- (b^2)/((sig^2)*(Q^2)))*(-(emu^2)/(2*(sig^2))) - (1/2)*log(sig^2 + (b^2)*(sigmax^2)) + log(pnorm(Q*(LOD_X - mux - (b*emu)/((sig^2)*(Q^2))))))

}

#' Log likelihood for a normally distributed covariate with a left censored outcome
#'
#' Calculates the log-likelihood contribution for a normally distriubted covariate when the outcome is left censored but the covariate is above the limit of detection
#'
#' @param X a normally distributed covariate
#' @param LOD_Y the limit of detection for the outcome
#' @param sig the standard deviation of the residuals
#' @param mux mean of the normally distributed covariate
#' @param sigmax standard deviation of the normally distributed covariate
#' @param alpha intercept: the average value of Y when all covariates = 0
#' @param b slope: the average change in Y for a 1 unit change in X
#' @param C optional vector of other covariates (not subject to LODs) included in the model
#' @param beta_C regression coefficients associated with C
#'
#' @return log-likelihood
#'
#' @examples
#' LL3_normal(LOD_Y=1, X = 2, sig = 0.5, mux = 1, sigmax = 1, alpha = 0, b = 1)
#'
#' @export
LL3_normal <- function(X, LOD_Y, mux, sigmax, sig, alpha, b, C=NULL, beta_C=NULL){

  if(is.null(C)==F){mu_y <- alpha + b*X + c(as.matrix(C)%*%beta_C)
  }else{mu_y <- alpha + b*X}


  return(dnorm(X,mean=mux, sd=sigmax, log=T) + pnorm(LOD_Y, mean = mu_y, sd=sig, log=T))

}

#' Log likelihood for a normally distributed covariate with a left censored outcome
#'
#' Calculates the log-likelihood contribution for a normally distriubted covariate when both the outcome and covariate are left censored.
#'
#' @param LOD_X the limit of detection for the normally distributed covariate
#' @param LOD_Y the limit of detection for the outcome
#' @param sig the standard deviation of the residuals
#' @param mux mean of the normally distributed covariate
#' @param sigmax standard deviation of the normally distributed covariate
#' @param alpha intercept: the average value of Y when all covariates = 0
#' @param b slope: the average change in Y for a 1 unit change in X
#' @param C optional vector of other covariates (not subject to LODs) included in the model
#' @param beta_C regression coefficients associated with C
#'
#' @return log-likelihood
#'
#' @examples
#' LL4_normal(LOD_Y=1, LOD_X = 0.2, sig = 0.5, mux = 1, sigmax = 1, alpha = 0, b = 1)
#'
#' @export
LL4_normal <- function(LOD_X, LOD_Y, mux, sigmax, sig, alpha, b, C=NULL, beta_C=NULL, single_LODs=F){

  if(single_LODs==T & is.null(C)==T){
    temp <- integrate(function(z){exp(LL2_normal(Y=z, LOD_X=unique(LOD_X), mux=mux, sigmax=sigmax, sig=sig, alpha=alpha, b=b, C=NULL, beta_C=NULL))}, lower = -Inf, upper = unique(LOD_Y))

    return(log(temp$value))
  }else{
    if(is.null(C)==F){
      temp <- sapply(seq(1, nrow(C),1), function(q){
        integrate(function(z){exp(LL2_normal(Y=z, LOD_X=LOD_X[q], mux=mux, sigmax=sigmax, sig=sig, alpha=alpha, b=b, C=(C[q,]), beta_C=beta_C))}, lower = -Inf, upper = LOD_Y[q])$value
      })

      return(log(temp))
    }else{

      temp <- sapply(seq(1, length(LOD_X),1), function(q){
        integrate(function(z){exp(LL2_normal(Y=z, LOD_X=LOD_X[q], mux=mux, sigmax=sigmax, sig=sig, alpha=alpha, b=b, C=NULL, beta_C=NULL))}, lower = -Inf, upper = LOD_Y[q])$value
      })

      return(log(temp))
    }
  }
}


#' Log likelihood for a log-normally distributed covariate without covariate or outcome censoring
#'
#' Calculates the log-likelihood contribution for a log-normally distriubted covariate when neither the covariate or outcome are censored.
#'
#' @param X a log-normally distributed covariate
#' @param Y outcome
#' @param sig the standard deviation of the residuals
#' @param mux mean of the ln(X)
#' @param sigmax standard deviation of ln(X)
#' @param alpha intercept: the average value of Y when all covariates = 0
#' @param b slope: the average change in Y for a 1 unit change in X
#' @param C optional vector of other covariates (not subject to LODs) included in the model
#' @param beta_C regression coefficients associated with C
#'
#' @return log-likelihood
#'
#' @examples
#' LL1_lognormal(X=2, Y=3, sig = 0.5, mux = 1, sigmax = 1, alpha = 0, b = 1)
#'
#' @export
LL1_lognormal <- function(X, Y, mux, sigmax, sig, alpha, b, beta_C=NULL, C=NULL){
  if(is.null(C)==F){dnorm(Y, mean = alpha+b*X+c(as.matrix(C)%*%beta_C), sd=sig, log=T)+dnorm(log(X), mean=mux, sd=sigmax, log=T)
  }else{dnorm(Y, mean = alpha+b*X, sd=sig, log=T)+dnorm(log(X), mean=mux, sd=sigmax, log=T)}
}

#' Log likelihood for a log-normally distributed left-censored covariate without outcome censoring
#'
#' Calculates the log-likelihood contribution for a log-normally distriubted covariate when the covariate is left censored and the outcome is not censored.
#'
#' @param LOD_X the limit of detection for the log-normally distributed covariate
#' @param Y outcome
#' @param sig the standard deviation of the residuals
#' @param mux mean of ln(X)
#' @param sigmax the standard deviation of ln(X)
#' @param alpha intercept: the average value of Y when all covariates = 0
#' @param b slope: the average change in Y for a 1 unit change in X
#' @param C optional vector of other covariates (not subject to LODs) included in the model
#' @param beta_C regression coefficients associated with C
#'
#' @return log-likelihood
#'
#' @examples
#' LL2_lognormal(Y=3, LOD_X = 0.2, sig = 0.5, mux = 1, sigmax = 1, alpha = 0, b = 1)
#'
#' @export
LL2_lognormal <- function(Y, LOD_X, mux, sigmax, sig, alpha, b, beta_C=NULL, C=NULL){
  if(is.null(C)==F){
    temp <- sapply(seq(1, length(LOD_X),1), function(q){
      integrate(function(z){exp(LL1_lognormal(X=z, Y=Y[q], mux=mux, sigmax=sigmax, sig=sig, alpha=alpha, b=b, beta_C=beta_C, C=C[q,]))}, lower = 0, upper = LOD_X[q])$value
    })
  }else{
    temp <- sapply(seq(1, length(LOD_X),1), function(q){
      integrate(function(z){exp(LL1_lognormal(X=z, Y=Y[q], mux=mux, sigmax=sigmax, sig=sig, alpha=alpha, b=b, beta_C=NULL, C=NULL))}, lower = 0, upper = LOD_X[q])$value
    })
  }

  return(log(temp))
}

#' Log likelihood for a log-normally distributed covariate without censoring when the outcome is below the limit of detection
#'
#' Calculates the log-likelihood contribution for a log-normally distriubted covariate when the covariate is uncensored and the outcome is left censored.
#'
#' @param X the log-normally distributed covariate
#' @param LOD_Y the limit of detection of the outcome
#' @param sig the standard deviation of the residuals
#' @param mux mean of ln(X)
#' @param sigmax standard deviation of ln(X)
#' @param alpha intercept: the average value of Y when all covariates = 0
#' @param b slope: the average change in Y for a 1 unit change in X
#' @param C optional vector of other covariates (not subject to LODs) included in the model
#' @param beta_C regression coefficients associated with C
#'
#' @return log-likelihood
#'
#' @examples
#' LL3_lognormal(LOD_Y=1, X = 2, sig = 0.5, mux = 1, sigmax = 1, alpha = 0, b = 1)
#'
#' @export
LL3_lognormal <- function(X, LOD_Y, mux, sigmax, sig, alpha, b, beta_C=NULL, C=NULL){
  if(is.null(C)==F){
    temp <- sapply(seq(1, length(X),1), function(q){
      integrate(function(z){exp(LL1_lognormal(X=X[q], Y=z, mux=mux, sigmax=sigmax, sig=sig, alpha=alpha, b=b, beta_C=beta_C, C=C[q,]))}, lower = -Inf, upper = LOD_Y[q])$value
    })
  }else{
    temp <- sapply(seq(1, length(X),1), function(q){
      integrate(function(z){exp(LL1_lognormal(X=X[q], Y=z, mux=mux, sigmax=sigmax, sig=sig, alpha=alpha, b=b, beta_C=NULL, C=NULL))}, lower = -Inf, upper = LOD_Y[q])$value
    })
  }
  return(log(temp))

}

#' Log likelihood for a log-normally distributed covariate when both the covariate and outcome are below the limit of detection
#'
#' Calculates the log-likelihood contribution for a log-normally distriubted covariate when the covariate and the outcome are both left censored.
#'
#' @param LOD_X the limit of detection of the log-normally distributed covariate
#' @param LOD_Y the limit of detection of the outcome
#' @param sig the standard deviation of the residuals
#' @param mux mean of ln(X)
#' @param sigmax standard deviation of ln(X)
#' @param alpha intercept: the average value of Y when all covariates = 0
#' @param b slope: the average change in Y for a 1 unit change in X
#' @param C optional vector of other covariates (not subject to LODs) included in the model
#' @param beta_C regression coefficients associated with C
#'
#' @return log-likelihood
#'
#' @examples
#' LL4_lognormal(LOD_Y=1, LOD_X = 0.2, sig = 0.5, mux = 1, sigmax = 1, alpha = 0, b = 1)
#'
#' @export
LL4_lognormal <- function(LOD_X, LOD_Y, mux, sigmax, sig, alpha, b, beta_C=NULL, C=NULL, single_LODs=F){

  if(is.null(C)==T & single_LODs==T){
    temp <- integrate(Vectorize(function(zz){
      integrate((function(z){
        dnorm(zz, mean = alpha+b*z, sd=sig, log=F)*dnorm(log(z), mean=mux, sd=sigmax, log=F)
      }), lower = 0, upper = unique(LOD_X))$value
    }), lower=-Inf, upper=unique(LOD_Y))$value

    return(log(temp))
  }else{
    if(is.null(C)==F){
      temp <- sapply(seq(1, nrow(C),1), function(q){
        integrate(Vectorize(function(zz){
          integrate((function(z){
            dnorm(zz, mean = alpha+b*z+c(as.matrix(C[q,])%*%beta_C), sd=sig, log=F)*dnorm(log(z), mean=mux, sd=sigmax, log=F)
          }), lower = 0, upper = LOD_X[1])$value
        }), lower=-Inf, upper=LOD_Y[q])$value
      })

      return(log(temp))
    }else{

      temp <- sapply(seq(1, length(LOD_X),1), function(q){
        integrate(Vectorize(function(zz){
          integrate((function(z){
            dnorm(zz, mean = alpha+b*z, sd=sig, log=F)*dnorm(log(z), mean=mux, sd=sigmax, log=F)
          }), lower = 0, upper = LOD_X[1])$value
        }), lower=-Inf, upper=LOD_Y[q])$value
      })

      return(log(temp))
    }
  }
}
