## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=F------------------------------------------------------------------
#  install.packages("devtools")

## ---- eval = F----------------------------------------------------------------
#  library(devtools)

## ---- eval = F----------------------------------------------------------------
#  install_github("camillemmoore/mle.lod)

## -----------------------------------------------------------------------------
library(mle.lod)
library(stats4)

## -----------------------------------------------------------------------------
set.seed(1000)
x <- exp(rnorm(100, mean=0, sd = 1))
y <- 10 + 2*x + rnorm(100, mean=0, sd=1)
plot(x,y)
x_censored <- ifelse(x < 0.5, 0.5, x)

## -----------------------------------------------------------------------------
summary(lm.1 <- lm(y ~ x))
m1.1 <- mle_fit(y = y, x = x_censored, 
              censor_x = ifelse(x_censored==0.5, 1, 0),
              censor_y = rep(0, length(y)), 
              x.dist = 'lognormal',
              start = c(mu_x = mean(log(x)), 
                        sigma_x = sd(log(x)), 
                        sigma_y = sd(lm.1$residuals), 
                        intercept = lm.1$coefficients[1], 
                        slope = lm.1$coefficients[2]))

summary(m1.1)

## -----------------------------------------------------------------------------
m1.2 <- mle_fit(y = y, x = x, 
              censor_x = ifelse(x < 0.5, 1, 0),
              LOD_x = 0.5,
              censor_y = rep(0, length(y)), 
              x.dist = 'lognormal',
              start = c(mu_x = mean(log(x)), 
                        sigma_x = sd(log(x)), 
                        sigma_y = sd(lm.1$residuals), 
                        intercept = lm.1$coefficients[1], 
                        slope = lm.1$coefficients[2]))

summary(m1.2)

## -----------------------------------------------------------------------------
x <- rnorm(100, mean=5, sd = 1)
y <- 10 + 2*x + rnorm(100, mean=0, sd=0.5)
plot(x,y)
LOD_x <- sample(c(3,4), size = 100, replace=T)
censor_x <- ifelse(x < LOD_x, 1, 0)

## -----------------------------------------------------------------------------
m2.1 <- mle_fit(y = y, x = x, 
              LOD_x = LOD_x,
              censor_x = censor_x,
              censor_y = rep(0, length(y)), 
              x.dist = 'normal')

summary(m2.1)

## -----------------------------------------------------------------------------
censor_y <- ifelse(y < 18, 1, 0)

m2.2 <- mle_fit(y = y, x = x, 
              LOD_x = LOD_x,
              LOD_y = 18,
              censor_x = censor_x,
              censor_y = censor_y, 
              x.dist = 'normal')

summary(m2.2)

## -----------------------------------------------------------------------------

summary(lm.2 <- lm(y~x))
m2.3 <- mle_fit(y = y, x = x, 
              LOD_x = LOD_x,
              LOD_y = 18,
              censor_x = censor_x,
              censor_y = censor_y, 
              x.dist = 'normal',
              start = c( mu_x = mean(x), 
                         sigma_x = sd(x), 
                         sigma_y = sd(lm.2$residuals), 
                         intercept = lm.2$coefficients[1], 
                         slope = lm.2$coefficients[2]))

summary(m2.3)

## -----------------------------------------------------------------------------

C <- data.frame(age = runif(100, 30, 50), gender = rbinom(size = 1, n=100, prob=0.5))
y <- 10 + 2*x + C$age*0.5 + C$gender*3 + rnorm(100, mean=0, sd=0.5)

m2.4 <- mle_fit(y = y, x = x, 
              LOD_x = LOD_x,
              censor_x = censor_x,
              censor_y = rep(0, length(y)), 
              covariates = C)

summary(m2.4)

