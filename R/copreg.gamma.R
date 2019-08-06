#' Bivariate Copula Regression with Gamma Marginal Distributions
#'
#' Bivariate copula regression models where the marginal distributions are of Gamma GLMs.
#'
#' @param f1 A regression formula for the first marginal GLM.
#' @param f2 A regression formula for the second marginal GLM.
#' @param copula A copula specified from the \code{copula} package.
#' @param data A matrix or data frame of observations. Categorical variables are allowed as covariates.
#' @param optim.method Method for the \code{optim} function.
#' @param optim.control Parameters that controls the \code{optim} function used for maximum likelihood.
#'
#' @return An object of class \code{copreg} providing the estimation results.
#'   The details of the output components are:
#'   \item{call}{The matched call.}
#'   \item{coefficients}{The estimated coefficients.}
#'   \item{copula.param}{The estimated copula parameters.}
#'   \item{copula}{The specified copula.}
#'   \item{fitted.values}{The fitted values of the regression.}
#'   \item{residuals}{The residuals from fitted values.}
#'   \item{loglike}{The final estimated maximum log-likelihood value.}
#'   \item{df}{Number of estimated parameters.}
#'   \item{AIC}{AIC values.}
#'   \item{BIC}{BIC values.}
#'   \item{formula}{The formulas used in the regression.}
#'   \item{y}{The input response data.}
#'   \item{n}{The number of observations in the data.}
#'   \item{Model.Matrix}{The used model matrix for each regression formula.}
#'
#' @examples
#' data("simdat.mcgr")
#' res <- copreg.gamma(f1 = y1 ~ x1+x2,
#'                     f2 = y2 ~ x1+x2,
#'                     copula = copula::gumbelCopula(dim=2),
#'                     data = simdat.mcgr)
#'
#' @export copreg.gamma

copreg.gamma <- function(f1,
                         f2,
                         copula,
                         data,
                         optim.method = "BFGS",
                         optim.control= list(fnscale=-1, trace=1)){
  tempcall<-as.call(c(expression(copreg.gamma),list(f1      = f1,
                                                    f2      = f2,
                                                    copula  = copula,
                                                    data    = substitute(data),
                                                    optim.method=optim.method,
                                                    optim.control=optim.control)))
  n <- nrow(data)
  namey1 <- as.character(f1[2])
  namey2 <- as.character(f2[2])
  y1     <- data[,names(data)==namey1]
  y2     <- data[,names(data)==namey2]
  y      <- cbind(y1,y2)
  newdata <- data
  newdata <- newdata[, names(newdata)!=namey1]
  newdata <- newdata[, names(newdata)!=namey2]
  if (!is.null(f2) && class(f1)=="formula" && as.character(f1[3])=="."){
    f1 <- stats::formula( paste( "y1", paste( names(newdata),'',collapse='+',sep='' ), sep='~'))
  }
  if (!is.null(f2) && class(f2)=="formula" && as.character(f2[3])=="."){
    f2 <- stats::formula( paste( "y2", paste( names(newdata),'',collapse='+',sep='' ), sep='~'))
  }
  xmat   <- list(stats::model.matrix(f1, data = data),
                 stats::model.matrix(f2, data = data))
  thecop <- copula
  #-------------------------------------------------------------
  # initial values
  m1 <- stats::glm(f1, family=Gamma(link="log"), data = data)
  m2 <- stats::glm(f2, family=Gamma(link="log"), data = data)
  fittedCop <- copula::fitCopula(thecop, copula::pobs(y))
  start.val <- as.vector(c(m1$coefficients,
                         log(MASS::gamma.shape(m1)$alpha),
                         m2$coefficients,
                         log(MASS::gamma.shape(m2)$alpha),
                         fittedCop@estimate))
  #-------------------------------------------------------------
  fit <- stats::optim(start.val, theloglik.gamma,
                      control = optim.control,
                      method = optim.method,
                      y = y, xmat = xmat, copula = thecop)
  finparam <- fit$par
  l1 <- ncol(xmat[[1]]) + 1
  l2 <- ncol(xmat[[2]]) + 1
  coef1 <- finparam[1:(l1-1)]
  coef2 <- finparam[(l1+1):(l1+l2-1)]
  pcop <- finparam[-(1:(l1 + l2))]
  shape1 <- exp(finparam[l1])
  shape2 <- exp(finparam[(l1+l2)])
  thecop@parameters <- pcop
  c1.fit <- exp( xmat[[1]] %*% coef1 )
  c2.fit <- exp( xmat[[2]] %*% coef2 )
  y1.residual <- y1 - c1.fit
  y2.residual <- y2 - c2.fit
  noparams <- length(finparam)
  AIC      <- -2*fit$value + noparams*2
  BIC      <- -2*fit$value + noparams*log(n)

  result <- list(coefficients   = c(beta1 =list(coef1),
                                    shape1=list(shape1),
                                    beta2 =list(coef2),
                                    shape2=list(shape2)),
                 copula.param   = pcop,
                 copula         = thecop,
                 fitted.values  = cbind(c1.fit, c2.fit),
                 AIC            = AIC,
                 BIC            = BIC,
                 loglike        = fit$value,
                 df             = noparams,
                 n              = n,
                 y              = y,
                 residuals      = cbind(y1.residual, y2.residual),
                 Model.Matrix   = xmat,
                 formula        = list(f1 = f1,
                                       f2 = f2),
                 call           = tempcall,
                 univariate.res = list(m1, m2))
  structure(result, class = 'copreg.gamma')
}



loglik.gamma.margin <- function(param, y, x) {
  l <- length(param)
  eta <- x %*% param[-l]
  sum(stats::dgamma(y, shape = exp(param[l]), rate = exp(param[l])/exp(eta), log = TRUE))
}
loglik.cop <- function(param, u, copula) {
  copula@parameters <- param
  sum( copula::dCopula(u, copula, log=TRUE ))
}
probtrans.gamma.margin <- function(param, y, x) {
  l <- length(param)
  eta <- x %*% param[-l]
  stats::pgamma(y, shape = exp(param[l]), rate = exp(param[l])/exp(eta))
}
theloglik.gamma <- function(theta, y, xmat, copula) {
  l1 <- ncol(xmat[[1]]) + 1
  l2 <- ncol(xmat[[2]]) + 1
  param1 <- theta[1:l1]
  param2 <- theta[(l1 + 1):(l1 + l2)]
  param.cop <- theta[-(1:(l1 + l2))]

  switch(class(copula),
         tCopula      = { cop.param.condition <- ifelse((param.cop[1]<=1 && param.cop[1]>=-1), TRUE, FALSE) },
         normalCopula = { cop.param.condition <- ifelse((param.cop<=1 && param.cop>=-1), TRUE, FALSE) },
         gumbelCopula = { cop.param.condition <- ifelse((param.cop>=1), TRUE, FALSE) },
         joeCopula    = { cop.param.condition <- ifelse((param.cop>=1), TRUE, FALSE) },
         claytonCopula= { cop.param.condition <- ifelse((param.cop>=-1 && param.cop!=0), TRUE, FALSE) },
         frankCopula  = { cop.param.condition <- ifelse((param.cop!=0), TRUE, FALSE) }
         )

  if ( cop.param.condition ){
    u <- cbind(probtrans.gamma.margin(param1, y[,1], xmat[[1]]),
               probtrans.gamma.margin(param2, y[,2], xmat[[2]]))
    copula@parameters <- param.cop
    loglik <- loglik.gamma.margin(param1, y[,1], xmat[[1]]) +
              loglik.gamma.margin(param2, y[,2], xmat[[2]]) +
              loglik.cop(param.cop, u, copula)
  } else {
    loglik <- NA
  }
  return(loglik)
}
