#' Maximum Likelihood Models with Censored Outcomes and Covariates
#'
#' Calculates model parameters when either the covariate, outcome or both are subject to a lower limit of detection.
#'
#' @param y a vector of potentially censored outcome values
#' @param x a vector of potentially censored covariate values
#' @param censor_x a vector the length of x; 0 if x is observed, 1 if censored
#' @param censor_y a vector the length of y; 0 if y is observed, 1 if censored
#' @param LOD_x a single LOD for the censored covariate or a vector of LODs (specific to each observation)
#' @param LOD_y a single LOD for the outcome or a vector of LODs (specific to each observation)
#' @param covariates an optional matrix or dataframe of additional uncensored covariates to include in the model
#' @param x.dist the distribution for censored predictor - either 'gaussian', 'normal' or 'lognormal'
#' @param start a vector of initial values for optimizer in the following order: mu_x, sigma_x, sigma_y, intercept, slope, covariates (in same order supplied in covariates.  If not supplied will default to mu_x = mean(x) if gaussian x, mean(log(x)) if lognormal x, sigma_x = sd(x) if gaussian x, sd(log(x)) if lognormal x, sigma_y = sd(y), intercept = mean(y), slope = 0, covariates = 0
#' @param method optimization method to use. See optim.
#' @param maxit integer. Maximum number of iterations for each optmization run.
#' @param reltol relative tolerance for stopping optimization.  See optim.
#' @param maxoptimit maximum number of opitimization iterations
#' @param llabstol absolute tolerance (change in log likelihood) for declaring convergence
#'
#' @return an object of class mle
#'
#' @examples
#' x <- rnorm(100, mean = 1)
#' y <- 3*x + rnorm(100, mean = 0, sd=0.2)
#' censor_x <- ifelse(x < 0.2, 1, 0)
#' censor_y <- ifelse(y < 0, 1, 0)
#' mle_fit(y=y, x=x,
#'        censor_x = censor_x, censor_y = censor_y,
#'        LOD_x = rep(0.2,100), LOD_y = rep(1,100),
#'        x.dist = 'normal')
#'
#' @export

mle_fit <- function(y, # censored outcome
                    x, # censored predictor
                    censor_x, # a vector the length of x; 0 if x is observed, 1 if censored
                    censor_y, # a vector the length of y; 0 if y is observed, 1 if censored
                    LOD_x=NULL, # a single LOD for the censored covariate or a vector of LODs (specific to each observation)
                    LOD_y=NULL, # a single LOD or a vector of LODs (specific to each observation)
                    covariates=NULL, # matrix or dataframe of covariates to include
                    x.dist = 'gaussian', # distribution for censored predictor - either gaussian or lognormal (will also accept normal),
                    start	= NULL, # Vector of initial values for optimizer in the following order:
                    #mu_x
                    #sigma_x
                    #sigma_y
                    #intercept
                    #slope
                    #covariates (in same order supplied in covariates
                    # if not supplied will default to
                    #mu_x = mean(x) if gaussian x, mean(log(x)) if lognormal x
                    #sigma_x = sd(x) if gaussian x, sd(log(x)) if lognormal x
                    #sigma_y = sd(y)
                    #intercept = mean(y)
                    #slope = 0
                    #covariates = 0
                    method = "Nelder-Mead", # Optimization method to use. See optim.
                    maxit =5000, # maximum number of iterations for each optmization run.
                    reltol = 1e-8, # relative tolerance for stopping optimization.  See optim.
                    maxoptimit = 10, # maximum number of opitimization iterations
                    llabstol = 1e-8 # absolute tolerance (change in log likelihood) for declaring convergence.
){
  #################################################################
  # CHECK INPUTS FOR ANY PROBLEMS AND DATA MANAGEMENT
  ################################################################
  # Check if we will use covariates
  use_covar <- !is.null(covariates)

  # Check for starting values
  if(is.null(start)==T){
    if(use_covar==F){start <- c(mu_x = ifelse(x.dist=='lognormal', mean(log(x)),mean(x)),
                                sigma_x = ifelse(x.dist=='lognormal', sd(log(x)),sd(x)),
                                sigma_y = sd(y),
                                intercept = mean(y),
                                slope = 0
    )
    }else{start <- c(mu_x = ifelse(x.dist=='lognormal', mean(log(x)),mean(x)),
                     sigma_x = ifelse(x.dist=='lognormal', sd(log(x)),sd(x)),
                     sigma_y = sd(y),
                     intercept = mean(y),
                     slope = 0, rep(0, ncol(covariates)))

    names(start)[6:(6+ncol(covariates)-1)] <- colnames(covariates)
    }
  }

  if((use_covar==F & length(start)!=5)){
    stop('The length of the starting values does not match the number of model paramters.')
  }

  if(use_covar==T){if(length(start)!=(5+ncol(covariates))){
    stop('The length of the starting values does not match the number of model paramters.')}
    if(nrow(covariates) != length(y)){stop('Number of covariate observations does not match number of y observations.')}
  }

  if(is.null(names(start))){stop('Provide names for starting values.')}

  # Check that dist is either gaussian or lognormal
  if((x.dist %in% c('gaussian', 'normal', 'lognormal'))==F){
    stop("Unavailable distribution.  x.dist must either be gaussian or lognormal.")}

  # Check for missing and that length of y and x are the same
  if(sum(is.na(y))>0){stop("There are missing values in y.")}
  if(sum(is.na(x))>0){stop("There are missing values in x.")}
  if(use_covar==T & anyNA(covariates)==TRUE){stop('There are missing covariate values.')}

  if(length(y) != length(x)){stop("x and y are not the same length.")}

  # Check that if censor is not null it is length of x and y, no missing values
  if(is.null(censor_y)==F){
    if(length(censor_y) != length(y)){stop("censor_y must be the same length as y.")}
    if(sum(is.na(censor_y))>0){stop("There are missing values in censor_y.")}
  }

  if(is.null(censor_x)==F){
    if(length(censor_x) != length(x)){stop("censor_x must be the same length as x.")}
    if(sum(is.na(censor_x))>0){stop("There are missing values in censor_x.")}
  }

  # If no censoring provided, assume all are observed
  if(is.null(censor_y)==T){censor_y <- rep(0, length(y))
  message('censor_y not provided.  Assuming all y values are observed.')
  }
  if(is.null(censor_x)==T){censor_x <- rep(0, length(x))
  message('censor_x not provided.  Assuming all x values are observed.')
  }

  # If no LOD provided, use observed value (assume that if marked censored, LOD value is provided)
  if(is.null(LOD_y)==T){LOD_y <- y
  message('No LOD for Y provided.  Assuming LOD is provided in y vector for censored values.')
  }
  if(is.null(LOD_x)==T){LOD_x <- x
  message('No LOD for X provided.  Assuming LOD is provided in x vector for censored values.')
  }

  # If single LOD provided, convert to a vector
  if(length(LOD_y)==1){LOD_y <- rep(LOD_y, length(y))
  message('Single LOD for Y provided.  Assuming all observations have the same limit of detection.')
  }
  if(length(LOD_x)==1){LOD_x <- rep(LOD_x, length(x))
  message('Single LOD for X provided.  Assuming all observations have the same limit of detection.')
  }

  # Determine if there is a single x and y LOD
  single_LODs <- (length(unique(LOD_x)) + length(unique(LOD_y))) == 2


  #################################################################
  # FIT THE MODEL
  ################################################################

  optim_it <- 1
  ll_diff <- ifelse(maxoptimit==1, NA, Inf)
  temp1 <- mle2(start = start, minuslogl = function(param){NNL(param, LOD_x=LOD_x, LOD_y=LOD_y, x=x, y=y, censor_x=censor_x, censor_y=censor_y, C=covariates, x.dist=x.dist)}, method = method, maxit = maxit, reltol = reltol)

  while(optim_it < maxoptimit & ll_diff > llabstol){
    temp2 <- mle2(start = coef(temp1),
                  minuslogl = function(param){NNL(param, LOD_x=LOD_x, LOD_y=LOD_y, x=x, y=y,
                                                  censor_x=censor_x, censor_y=censor_y, C=covariates,
                                                  x.dist=x.dist)}, method = method, maxit = maxit, reltol = reltol)

    ll_diff <- logLik(temp2)-logLik(temp1)
    optim_it <- optim_it + 1
    temp1 <- temp2

  }

  if(maxoptimit==1){message('Only a single optimization iteration allowed.  Convergence not guaranteed.')}
  if(maxoptimit!=1 & ll_diff > llabstol){warning('Model did not converge in allowed number of optimization interations.')}

  attributes(temp1)$ll_diff <- ll_diff
  attributes(temp1)$optim_it <- optim_it
  temp1

}
