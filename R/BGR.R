#' Bivariate Gamma Regression
#'
#' Bivariate gamma regression models for three model types: EE, EI, IE. Models are estimated by EM algorithm.
#'
#' @param modelName A character string indicating which model to be fitted. Need to be one of "EE", "EI", "IE".
#' @param y A vector of character strings indicating which variables in the data are treated as response or dependent variables.
#' @param data A matrix or data frame of observations. Categorical variables are allowed as covariates.
#' @param f1 A regression formula for the \eqn{latex}{\alpha_1} parameter in the bivariate gamma distribution. Note that, depending on the model type, might not be necessary to provide it.
#' @param f2 A regression formula for the \eqn{latex}{\alpha_2} parameter in the bivariate gamma distribution. Note that, depending on the model type, might not be necessary to provide it.
#' @param f3 A regression formula for the \eqn{latex}{\alpha_3} parameter in the bivariate gamma distribution. Note that, depending on the model type, might not be necessary to provide it.
#' @param f4 A regression formula for the \eqn{latex}{\beta} parameter in the bivariate gamma distribution. Note that, depending on the model type, might not be necessary to provide it.
#' @param maxit A parameter that controls the number of maximum iteration in the EM algorithm. The default is 100.
#' @param tol A parameter that controls the convergence tolerance in the EM algorithm. The default is 1e-5.
#' @param verbose A logical controlling whether estimations in each EM iteration are shown in the fitting process. The default is TRUE.
#'
#' @return An object of class \code{BGR} providing the estimation results.
#'   The details of the output components are:
#'   \item{call}{The matched call.}
#'   \item{coefficients}{The estimated coefficients.}
#'   \item{alpha1}{The estimated alpha1 values.}
#'   \item{alpha2}{The estimated alpha2 values.}
#'   \item{alpha3}{The estimated alpha3 values.}
#'   \item{beta}{The estimated beta values.}
#'   \item{fitted.values}{The fitted values of the regression.}
#'   \item{loglike}{The final estimated maximum log-likelihood value.}
#'   \item{ll}{The sequence of log-likelihood values in the EM algorithm fitting process.}
#'   \item{df}{Number of estimated parameters.}
#'   \item{AIC}{AIC values.}
#'   \item{BIC}{BIC values.}
#'   \item{iter}{Total iteration numbers.}
#'   \item{formula}{The formulas used in the regression.}
#'   \item{y}{The input response data.}
#'   \item{n}{The number of observations in the data.}
#'   \item{Model.Matrix}{The used model matrix for each regression formula.}
#'   \item{trace}{All estimated coefficients and alpha, beta values in the EM algorithm.}
#'
#' @examples
#' \donttest{
#' mod1 <- BGR(modelName = "EE",
#'             y = c("y1","y2"), data = fullsim,
#'             f1 = ~ w1 + w2,
#'             f2 = ~ w2 + w3,
#'             f3 = ~ w1 + w2 + w3,
#'             f4 = ~ w1 + w2 + w3,
#'             verbose= FALSE)
#' mod1
#'
#' mod2 <- BGR(modelName = "EI",
#'             y = c("y1","y2"), data = fullsim,
#'             f1     = ~ w1 + w2,
#'             f2     = ~ w2 + w3,
#'             f3     = ~ w1 + w2 + w3,
#'             verbose= FALSE)
#' mod2
#' mod3 <- BGR(modelName = "IE",
#'             y = c("y1","y2"), data = fullsim,
#'             f4     = ~ w1 + w2 + w3,
#'             verbose= FALSE)
#' mod3
#' }
#'
#' @importFrom stats glm optim uniroot model.matrix Gamma formula runif coef
#' @export

BGR <- function(modelName = c("EE","EI","IE"),
                y,
                data,
                f1,
                f2,
                f3,
                f4,
                maxit=100,
                tol=1e-5,
                verbose=TRUE){
  switch(modelName,
         EE = {    # BGR.EE: all alpha's and beta are regressed on covariates
           if (is.null(f1)){ stop("f1 must be supplied for EE model type") }
           if (is.null(f2)){ stop("f2 must be supplied for EE model type") }
           if (is.null(f3)){ stop("f3 must be supplied for EE model type") }
           if (is.null(f4)){ stop("f4 must be supplied for EE model type") }
           return(BGR_EE(data=data, y=y,
                         f1=f1, f2=f2, f3=f3, f4=f4,
                         maxit=maxit, tol=tol,
                         verbose=verbose))
           },
         EI = {    # BGR.EI: same beta for all observations
           if (is.null(f1)){ stop("f1 must be supplied for EI model type") }
           if (is.null(f2)){ stop("f2 must be supplied for EI model type") }
           if (is.null(f3)){ stop("f3 must be supplied for EI model type") }
           return(BGR_EI(data=data, y=y,
                         f1=f1, f2=f2, f3=f3,
                         maxit=maxit, tol=tol,
                         verbose=verbose))
           },
         IE = {# BGR.IE: only beta is regressed on its covariates, not alpha
           if (is.null(f4)){ stop("f4 must be supplied for IE model type") }
           return(BGR_IE(data=data, y=y,
                         f4=f4,
                         maxit=maxit, tol=tol,
                         verbose=verbose))
           },
         stop("invalid model type name")
         )
}

#' @rdname BGR
#' @export
BGR_EE <- function(y,
                   data,
                   f1,
                   f2,
                   f3,
                   f4,
                   maxit  = 100,
                   tol    = 1e-5,
                   verbose= TRUE){
  options(warn=-1)
  tempcall<-as.call(c(expression(BGR),list(y       = y,
                                           data    = substitute(data),
                                           f1      = f1,
                                           f2      = f2,
                                           f3      = f3,
                                           f4      = f4,
                                           maxit   = maxit,
                                           tol     = tol,
                                           verbose = verbose)))
  namey1 <- y[1]
  namey2 <- y[2]
  y1     <- data[,names(data)==namey1]
  y2     <- data[,names(data)==namey2]
  n      <- length(y1)

  newdata <- data;
  newdata <- newdata[, names(newdata)!=namey1];
  newdata <- newdata[, names(newdata)!=namey2]

  if (!is.null(f1)){
    if (as.character(f1[2])=="."){
      f1.new <- formula( paste( "x1", paste( names(newdata),'',collapse='+',sep='' ), sep='~'))
    } else {f1.new <- formula(paste("x1", as.character(f1)[2], sep="~"))}
  } else stop("f1 cannot be NULL, must be supplied")
  if (!is.null(f2)){
    if (as.character(f2[2])=="."){f2.new <- formula( paste( "x2", paste( names(newdata),'',collapse='+',sep='' ), sep='~'))}
    else {f2.new <- formula(paste("x2", as.character(f2)[2], sep="~"))}
  } else stop("f2 cannot be NULL, must be supplied")
  if (!is.null(f3)){
    if (as.character(f3[2])=="."){f3.new <- formula( paste( "x3", paste( names(newdata), "", collapse = "+", sep = ""), sep="~"))}
    else {f3.new <- formula(paste("x3", as.character(f3)[2], sep="~"))}
  } else stop("f3 cannot be NULL, must be supplied")
  if (!is.null(f4)){
    if (as.character(f4[2])=="."){f4.new <- formula( paste( "be", paste( names(newdata), "", collapse = "+", sep = ""), sep="~"))}
    else {f4.new <- formula(paste("be", as.character(f4)[2], sep="~"))}
  } else stop("f4 cannot be NULL, must be supplied")

  Model.Matrix.1   <- stats::model.matrix(f1.new[-2], data=newdata)
  Model.Matrix.2   <- stats::model.matrix(f2.new[-2], data=newdata)
  Model.Matrix.3   <- stats::model.matrix(f3.new[-2], data=newdata)
  Model.Matrix.4   <- stats::model.matrix(f4.new[-2], data=newdata)

  #-----------------------------
  # starting values
  x3.start <- rep(0, n)
  for (i in seq_len(n)){
    x3.start[i] <- runif(1, min=0, max=min(y1[i], y2[i]))
  }
  x1.start         <- y1 - x3.start
  x2.start         <- y2 - x3.start
  start.temp.glm.1 <- stats::glm(x1.start~Model.Matrix.1[,-1], family=Gamma(link="log"))
  start.temp.glm.2 <- stats::glm(x2.start~Model.Matrix.2[,-1], family=Gamma(link="log"))
  start.temp.glm.3 <- stats::glm(x3.start~Model.Matrix.3[,-1], family=Gamma(link="log"))
  coef1.current    <- start.temp.glm.1$coefficient
  coef2.current    <- start.temp.glm.2$coefficient
  coef3.current    <- start.temp.glm.3$coefficient
  alpha1.current   <- rep(1/summary(start.temp.glm.1)$dispersion, n)
  alpha2.current   <- rep(1/summary(start.temp.glm.2)$dispersion, n)
  alpha3.current   <- rep(1/summary(start.temp.glm.3)$dispersion, n)
  beta.current.1   <- 1/summary(start.temp.glm.1)$dispersion / stats::predict.glm(start.temp.glm.1, type="response")
  beta.current.2   <- 1/summary(start.temp.glm.2)$dispersion / stats::predict.glm(start.temp.glm.2, type="response")
  beta.current.3   <- 1/summary(start.temp.glm.3)$dispersion / stats::predict.glm(start.temp.glm.3, type="response")
  beta.current     <- rowMeans(cbind(beta.current.1, beta.current.2, beta.current.3))
  start.temp.glm.4 <- stats::glm(beta.current~Model.Matrix.4[,-1], family=Gamma(link="log"))
  coef4.current    <- start.temp.glm.4$coefficient

  #--------------------------
  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)
  j            <- 1

  while ( (loglike.diff > tol) && (j <= maxit) ) {

    #-------
    #E step
    #-------
    Expected.x3    <- rep(0,n)
    Expected.x1    <- rep(0,n)
    Expected.x2    <- rep(0,n)
    Expected.logx1 <- rep(0,n)
    Expected.logx2 <- rep(0,n)
    Expected.logx3 <- rep(0,n)
    den            <- rep(0,n)
    for (i in seq_len(n)){
      expected.latent.res <- expected_latent(c(y1[i],y2[i]),
                                             alpha=c(alpha1.current[i],alpha2.current[i],alpha3.current[i]),
                                             beta =beta.current[i])
      Expected.x3[i]      <- expected.latent.res$Expected.s
      Expected.x1[i]      <- y1[i]-expected.latent.res$Expected.s
      Expected.x2[i]      <- y2[i]-expected.latent.res$Expected.s
      Expected.logx3[i]   <- expected.latent.res$Expected.logs
      Expected.logx1[i]   <- expected.latent.res$Expected.logy1s
      Expected.logx2[i]   <- expected.latent.res$Expected.logy2s
      if (Expected.x1[i] <=0 || Expected.x2[i] <=0){
        Expected.x3[i]  <- min(y1[i], y2[i])/2
        Expected.x1[i]  <- y1[i]-Expected.x3[i]
        Expected.x2[i]  <- y2[i]-Expected.x3[i]
        warning("negative expected x1 or x2 at ", i, " obs at ", j, " iteration", "\n")
      }
      den[i]              <- expected.latent.res$denominator
      #if (expected.latent.res$numerator.integration.message != "OK"){print(expected.latent.res$numerator.integration.message)}
      #if (expected.latent.res$denominator.integration.message != "OK"){print(expected.latent.res$denominator.integration.message)}
    }

    loglike.new    <- sum(den)
    loglike[j]     <- loglike.new
    loglike.diff <- abs(loglike.new-loglike.current)/(1+abs(loglike.new))
    if (verbose){
      cat(c("iter:", j, "; current loglike:", round(loglike.new,4), "; loglike change:", round(loglike.diff,6), "\n"))
    }

    #--------
    # M step
    #--------
    beta.temp      <- (alpha1.current + alpha2.current + alpha3.current) / (Expected.x1 + Expected.x2 + Expected.x3)
    m4             <- stats::glm(beta.temp ~ Model.Matrix.4[,-1], family=Gamma(link="log"))
    Q.b.function   <- function(coef){
      q.b.res <- sum((alpha1.current + alpha2.current + alpha3.current) * log(exp(coef %*% t(Model.Matrix.4)))) -
        sum(exp(coef %*% t(Model.Matrix.4)) * (Expected.x1 + Expected.x2 + Expected.x3))
      return(q.b.res)
    }
    coef4.optim    <- stats::optim(par     = as.vector(m4$coefficients),
                            fn      = Q.b.function,
                            gr      = NULL,
                            hessian = TRUE,
                            control = list(fnscale=-1, maxit=1e+5))
    if (coef4.optim$convergence!=0) cat("optim coef. beta not converged")
    coef4.new      <- coef4.optim$par
    beta.new       <- as.vector(exp(coef4.new%*%t(Model.Matrix.4)))

    m1             <- stats::glm(alpha1.current~Model.Matrix.1[,-1], family=Gamma(link="log"))
    Q.function.a1  <- function(coef){
      q.res <- sum( exp(coef %*% t(Model.Matrix.1)) * log(beta.new) ) -
        sum( lgamma(exp(coef%*%t(Model.Matrix.1))) ) +
        sum( exp(coef %*% t(Model.Matrix.1)) * Expected.logx1 )
      return(q.res)
    }
    coef1.optim    <- stats::optim(par     = as.vector(m1$coefficients),
                            fn      = Q.function.a1,
                            gr      = NULL,
                            hessian = TRUE,
                            control = list(fnscale=-1, maxit=1e+5))
    if (coef1.optim$convergence!=0) cat("optim coef a1 not converged", "\n")
    coef1.new      <- coef1.optim$par
    alpha1.new     <- as.vector(exp(coef1.new%*%t(Model.Matrix.1)))

    Q.function.a2  <- function(coef){
      q.res <- sum( exp(coef %*% t(Model.Matrix.2)) * log(beta.new) ) -
        sum( lgamma(exp(coef%*%t(Model.Matrix.2))) ) +
        sum( exp(coef %*% t(Model.Matrix.2)) * Expected.logx2 )
      return(q.res)
    }
    m2             <- stats::glm(alpha2.current~Model.Matrix.2[,-1], family=Gamma(link="log"))
    coef2.optim    <- stats::optim(par     = as.vector(m2$coefficients),
                            fn      = Q.function.a2,
                            gr      = NULL,
                            hessian = TRUE,
                            control = list(fnscale=-1, maxit=1e+5))
    if (coef2.optim$convergence) cat("optim coef a2 not converged", "\n")
    coef2.new      <- coef2.optim$par
    alpha2.new     <- as.vector(exp(coef2.new%*%t(Model.Matrix.2)))

    Q.function.a3 <- function(coef){
      q.res <- sum( exp(coef %*% t(Model.Matrix.3)) * log(beta.new) ) -
        sum( lgamma(exp(coef%*%t(Model.Matrix.3))) ) +
        sum( exp(coef %*% t(Model.Matrix.3)) * Expected.logx3 )
      return(q.res)
    }
    m3             <- stats::glm(alpha3.current~Model.Matrix.3[,-1], family=Gamma(link="log"))
    coef3.optim    <- stats::optim(par     = as.vector(m3$coefficients),
                            fn      = Q.function.a3,
                            gr      = NULL,
                            hessian = TRUE,
                            control = list(fnscale=-1, maxit=1e+5))
    if (coef3.optim$convergence!=0) cat("optim coef a3 not converged", "\n")
    coef3.new <- coef3.optim$par
    alpha3.new     <- as.vector(exp(coef3.new%*%t(Model.Matrix.3)))

    if (verbose){
      cat(c("alpha1:", round(summary(alpha1.new),4), "\n"))
      cat(c("alpha2:", round(summary(alpha2.new),4), "\n"))
      cat(c("alpha3:", round(summary(alpha3.new),4), "\n"))
      cat(c("beta:",   round(summary(beta.new),4),  "\n"))
    }

    coef1.current  <- coef1.new
    coef2.current  <- coef2.new
    coef3.current  <- coef3.new
    coef4.current  <- coef4.new
    loglike.current<- loglike.new
    alpha1.current <- alpha1.new
    alpha2.current <- alpha2.new
    alpha3.current <- alpha3.new
    beta.current   <- beta.new
    j              <- j+1
  }

  if (!all(loglike[1:(j-1)] == cummax(loglike[1:(j-1)]))){
    cat("Loglikelihood is not monotonically increasing!", "\n")
  }

  noparams<- m1$rank + m2$rank + m3$rank + m3$rank
  AIC     <- -2*loglike[j-1] + noparams * 2
  BIC     <- -2*loglike[j-1] + noparams * log(n)

  # fitted values
  y1.fitted <- (alpha1.current + alpha3.current) / beta.current
  y2.fitted <- (alpha2.current + alpha3.current) / beta.current
  y1.residual <- y1 - y1.fitted
  y2.residual <- y2 - y2.fitted

  # formula
  f1.formula <- formula(paste("", as.character(f1.new)[3], sep="~"))
  f2.formula <- formula(paste("", as.character(f2.new)[3], sep="~"))
  f3.formula <- formula(paste("", as.character(f3.new)[3], sep="~"))
  f4.formula <- formula(paste("", as.character(f4.new)[3], sep="~"))

  result<-list(modelName   = "EE",
               coefficients= list(f1.coef=coef1.new,
                                  f2.coef=coef2.new,
                                  f3.coef=coef3.new,
                                  f4.coef=coef4.new),
               alpha1      = alpha1.current,
               alpha2      = alpha2.current,
               alpha3      = alpha3.current,
               beta        = beta.current,
               fitted.values= data.frame(y1 = y1.fitted,
                                         y2 = y2.fitted),
               residuals   = cbind(y1.residual, y2.residual),
               loglike     = loglike[j-1],
               ll          = loglike[1:(j-1)],
               df          = noparams,
               AIC         = AIC,
               BIC         = BIC,
               y           = cbind(y1, y2),
               n           = n,
               Model.Matrix= list(Model.Matrix.1,
                                  Model.Matrix.2,
                                  Model.Matrix.3,
                                  Model.Matrix.4),
               formula     = list(f1.formula,
                                  f2.formula,
                                  f3.formula,
                                  f4.formula),
               call        = tempcall,
               iter        = (j-1))
  options(warn=0)
  structure(result, class = c('BGR'))
}

#' @rdname BGR
#' @export
BGR_EI <- function(y,
                   data,
                   f1,
                   f2,
                   f3,
                   maxit  = 100,
                   tol    = 1e-5,
                   verbose= TRUE){
  options(warn=-1)
  tempcall <- as.call( c(expression(BGR), list(y       = y,
                                               data    = substitute(data),
                                               f1      = f1,
                                               f2      = f2,
                                               f3      = f3,
                                               maxit   = maxit,
                                               tol     = tol,
                                               verbose = verbose) ) )
  namey1 <- y[1]
  namey2 <- y[2]
  y1     <- data[,names(data)==namey1]
  y2     <- data[,names(data)==namey2]
  n      <- length(y1)

  newdata <- data
  newdata <- newdata[, names(newdata)!=namey1]
  newdata <- newdata[, names(newdata)!=namey2]

  if (!is.null(f1)){
    if (as.character(f1[2])=="."){
      f1.new <- formula( paste( "x1", paste( names(newdata),'',collapse='+',sep='' ), sep='~'))
    } else {f1.new <- formula(paste("x1", as.character(f1)[2], sep="~"))}
  } else stop("f1 cannot be NULL, must be supplied")
  if (!is.null(f2)){
    if (as.character(f2[2])=="."){f2.new <- formula( paste( "x2", paste( names(newdata),'',collapse='+',sep='' ), sep='~'))}
    else {f2.new <- formula(paste("x2", as.character(f2)[2], sep="~"))}
  } else stop("f2 cannot be NULL, must be supplied")
  if (!is.null(f3)){
    if (as.character(f3[2])=="."){f3.new <- formula( paste( "x3", paste( names(newdata), "", collapse = "+", sep = ""), sep="~"))}
    else {f3.new <- formula(paste("x3", as.character(f3)[2], sep="~"))}
  } else stop("f3 cannot be NULL, must be supplied")

  Model.Matrix.1 <- stats::model.matrix(f1.new[-2], data=newdata)
  Model.Matrix.2 <- stats::model.matrix(f2.new[-2], data=newdata)
  Model.Matrix.3 <- stats::model.matrix(f3.new[-2], data=newdata)
  #-----------------------------
  # starting values
  x3.start <- rep(0, n)
  for (i in seq_len(n)){
    x3.start[i] <- runif(1, min=0, max=min(y1[i], y2[i]))
  }
  x1.start       <- y1 - x3.start
  x2.start       <- y2 - x3.start
  start.temp.glm.1 <- stats::glm(x1.start~Model.Matrix.1[,-1], family=Gamma(link="log"))
  start.temp.glm.2 <- stats::glm(x2.start~Model.Matrix.2[,-1], family=Gamma(link="log"))
  start.temp.glm.3 <- stats::glm(x3.start~Model.Matrix.3[,-1], family=Gamma(link="log"))
  coef1.current  <- start.temp.glm.1$coefficient
  coef2.current  <- start.temp.glm.2$coefficient
  coef3.current  <- start.temp.glm.3$coefficient
  alpha1.current <- rep(1/summary(start.temp.glm.1)$dispersion, n)
  alpha2.current <- rep(1/summary(start.temp.glm.2)$dispersion, n)
  alpha3.current <- rep(1/summary(start.temp.glm.3)$dispersion, n)
  beta.current.1 <- 1/summary(start.temp.glm.1)$dispersion / stats::predict.glm(start.temp.glm.1, type="response")
  beta.current.2 <- 1/summary(start.temp.glm.2)$dispersion / stats::predict.glm(start.temp.glm.2, type="response")
  beta.current.3 <- 1/summary(start.temp.glm.3)$dispersion / stats::predict.glm(start.temp.glm.3, type="response")
  beta.current   <- mean(cbind(beta.current.1, beta.current.2, beta.current.3))
  #--------------------------

  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)
  j            <- 1

  while ( (loglike.diff > tol) && (j <= maxit) ) {

    #-------
    #E step
    #-------
    Expected.x3    <- rep(0,n)
    Expected.x1    <- rep(0,n)
    Expected.x2    <- rep(0,n)
    Expected.logx1 <- rep(0,n)
    Expected.logx2 <- rep(0,n)
    Expected.logx3 <- rep(0,n)
    den            <- rep(0,n)
    for (i in seq_len(n)){
      expected.latent.res <- expected_latent(c(y1[i],y2[i]),
                                             alpha=c(alpha1.current[i],alpha2.current[i],alpha3.current[i]),
                                             beta =beta.current)
      Expected.x3[i]      <- expected.latent.res$Expected.s
      Expected.x1[i]      <- y1[i]-expected.latent.res$Expected.s
      Expected.x2[i]      <- y2[i]-expected.latent.res$Expected.s
      Expected.logx3[i]   <- expected.latent.res$Expected.logs
      Expected.logx1[i]   <- expected.latent.res$Expected.logy1s
      Expected.logx2[i]   <- expected.latent.res$Expected.logy2s
      if (Expected.x1[i] <=0 || Expected.x2[i] <=0){
        Expected.x3[i]  <- min(y1[i], y2[i])/2
        Expected.x1[i]  <- y1[i]-Expected.x3[i]
        Expected.x2[i]  <- y2[i]-Expected.x3[i]
        warning("negative expected x1 or x2 at ", i, " obs at ", j, " iteration", "\n")
      }
      den[i]              <- expected.latent.res$denominator
      #if (expected.latent.res$numerator.integration.message != "OK"){print(expected.latent.res$numerator.integration.message)}
      #if (expected.latent.res$denominator.integration.message != "OK"){print(expected.latent.res$denominator.integration.message)}
    }

    loglike.new    <- sum(den)
    loglike[j]     <- loglike.new
    loglike.diff <- abs(loglike.new-loglike.current)/(1+abs(loglike.new))
    if (verbose){
      cat(c("iter:", j, "; current loglike:", round(loglike.new,4), "; loglike change:", round(loglike.diff,6), "\n"))
    }

    #--------
    # M step
    #--------

    Q.function <- function(coef, Model.Matrix, Expected.logs){
      q.res <- sum( exp(coef%*%t(Model.Matrix)) * log(beta.current) ) -
        sum( lgamma(exp(coef%*%t(Model.Matrix))) ) +
        sum( exp(coef%*%t(Model.Matrix)) * Expected.logs )
      return(q.res)
    }
    m1             <- stats::glm(alpha1.current~Model.Matrix.1[,-1], family=Gamma(link="log"))
    coef1.new      <- stats::optim(par          = as.vector(m1$coefficients),
                            fn           = Q.function,
                            Model.Matrix = Model.Matrix.1,
                            Expected.logs= Expected.logx1,
                            gr           = NULL,
                            hessian      = TRUE,
                            control      = list(fnscale=-1))$par
    alpha1.new     <- as.vector(exp(coef1.new%*%t(Model.Matrix.1)))

    m2             <- stats::glm(alpha2.current~Model.Matrix.2[,-1], family=Gamma(link="log"))
    coef2.new      <- stats::optim(par          = as.vector(m2$coefficients),
                            fn           = Q.function,
                            Model.Matrix = Model.Matrix.2,
                            Expected.logs= Expected.logx2,
                            gr           = NULL,
                            hessian      = TRUE,
                            control      = list(fnscale=-1))$par
    alpha2.new     <- as.vector(exp(coef2.new%*%t(Model.Matrix.2)))

    m3             <- stats::glm(alpha3.current~Model.Matrix.3[,-1], family=Gamma(link="log"))
    coef3.new      <- stats::optim(par          = as.vector(m3$coefficients),
                            fn           = Q.function,
                            Model.Matrix = Model.Matrix.3,
                            Expected.logs= Expected.logx3,
                            gr           = NULL,
                            hessian      = TRUE,
                            control      = list(fnscale=-1))$par
    alpha3.new     <- as.vector(exp(coef3.new%*%t(Model.Matrix.3)))

    beta.new       <- sum(alpha1.new + alpha2.new + alpha3.new) / sum(Expected.x1 + Expected.x2 + Expected.x3)

    if (verbose){
      cat(c("alpha1:", round(summary(alpha1.new),4), "\n"))
      cat(c("alpha2:", round(summary(alpha2.new),4), "\n"))
      cat(c("alpha3:", round(summary(alpha3.new),4), "\n"))
      cat(c("beta:",   round(beta.new,4) , "\n"))
    }

    loglike.current<- loglike.new
    coef1.current  <- coef1.new
    coef2.current  <- coef2.new
    coef3.current  <- coef3.new
    alpha1.current <- alpha1.new
    alpha2.current <- alpha2.new
    alpha3.current <- alpha3.new
    beta.current   <- beta.new
    j              <- j+1
  }

  if (!all(loglike[1:(j-1)] == cummax(loglike[1:(j-1)]))){
    cat("Loglikelihood is not monotonically increasing!", "\n")
  }
  noparams <- m1$rank + m2$rank + m3$rank + 1
  AIC      <- -2*loglike[j-1] + noparams * 2
  BIC      <- -2*loglike[j-1] + noparams * log(n)

  # fitted values
  y1.fitted <- (alpha1.current + alpha3.current) / beta.current
  y2.fitted <- (alpha2.current + alpha3.current) / beta.current
  y1.residual <- y1 - y1.fitted
  y2.residual <- y2 - y2.fitted

  # formula
  f1.formula <- formula(paste("", as.character(f1.new)[3], sep="~"))
  f2.formula <- formula(paste("", as.character(f2.new)[3], sep="~"))
  f3.formula <- formula(paste("", as.character(f3.new)[3], sep="~"))

  result<-list(modelName   = "EI",
               coefficients= list(f1.coef=coef1.new,
                                  f2.coef=coef2.new,
                                  f3.coef=coef3.new),
               alpha1      = alpha1.current,
               alpha2      = alpha2.current,
               alpha3      = alpha3.current,
               beta        = beta.current,
               fitted.values=data.frame(y1=y1.fitted,
                                        y2=y2.fitted),
               residuals   = cbind(y1.residual, y2.residual),
               loglike     = loglike[j-1],
               ll          = loglike[1:(j-1)],
               df          = noparams,
               AIC         = AIC,
               BIC         = BIC,
               y           = cbind(y1, y2),
               n           = n,
               Model.Matrix= list(Model.Matrix.1,
                                  Model.Matrix.2,
                                  Model.Matrix.3),
               formula     = list(f1.formula,
                                  f2.formula,
                                  f3.formula),
               call        = tempcall,
               iter        = j-1)
  options(warn=0)
  structure(result, class = c('BGR'))
}

#' @rdname BGR
#' @export
BGR_IE <- function(y,
                   data,
                   f4,
                   maxit  = 100,
                   tol    = 1e-5,
                   verbose= TRUE){
  options(warn=-1)
  tempcall<-as.call(c(expression(BGR),list(y       = y,
                                           data    = substitute(data),
                                           f4      = f4,
                                           maxit   = maxit,
                                           tol     = tol,
                                           verbose = verbose) ) )
  namey1 <- y[1]
  namey2 <- y[2]
  y1     <- data[,names(data)==namey1]
  y2     <- data[,names(data)==namey2]
  n      <- length(y1)

  newdata <- data
  newdata <- newdata[,names(newdata)!=namey1]
  newdata <- newdata[,names(newdata)!=namey2]

  if (!is.null(f4)){
    if (as.character(f4[2])=="."){f4.new <- formula( paste( "be", paste( names(newdata), "", collapse = "+", sep = ""), sep="~"))}
    else {f4.new <- formula(paste("be", as.character(f4)[2], sep="~"))}
  } else stop("f4 cannot be NULL, must be supplied")

  Model.Matrix.4   <- stats::model.matrix(f4.new[-2], data=newdata)
  #-----------------------------
  # starting values
  x3.start <- rep(0, n)
  for (i in seq_len(n)){
    x3.start[i] <- runif(1, min=0, max=min(y1[i], y2[i]))
  }
  x1.start         <- y1 - x3.start
  x2.start         <- y2 - x3.start

  alpha1.current   <- MASS::fitdistr(x1.start, "gamma")$estimate[1]
  alpha2.current   <- MASS::fitdistr(x1.start, "gamma")$estimate[1]
  alpha3.current   <- MASS::fitdistr(x1.start, "gamma")$estimate[1]
  beta.current     <- rowMeans(cbind(alpha1.current/x1.start, alpha2.current/x2.start, alpha3.current/x3.start))
  start.temp.glm   <- stats::glm(beta.current~Model.Matrix.4[,-1], family=Gamma(link="log"))
  coef4.current    <- start.temp.glm$coefficient
  #--------------------------

  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)
  j            <- 1

  while ( (loglike.diff > tol) && (j <= maxit) ) {
    #-------
    #E step
    #-------
    Expected.x3    <- rep(0,n)
    Expected.x1    <- rep(0,n)
    Expected.x2    <- rep(0,n)
    Expected.logx1 <- rep(0,n)
    Expected.logx2 <- rep(0,n)
    Expected.logx3 <- rep(0,n)
    den            <- rep(0,n)
    for (i in seq_len(n)){
      expected.latent.res <- expected_latent(c(y1[i],y2[i]),
                                             alpha=c(alpha1.current,alpha2.current,alpha3.current),
                                             beta =beta.current[i])
      Expected.x3[i]      <- expected.latent.res$Expected.s
      Expected.x1[i]      <- y1[i]-expected.latent.res$Expected.s
      Expected.x2[i]      <- y2[i]-expected.latent.res$Expected.s
      Expected.logx3[i]   <- expected.latent.res$Expected.logs
      Expected.logx1[i]   <- expected.latent.res$Expected.logy1s
      Expected.logx2[i]   <- expected.latent.res$Expected.logy2s
      if (Expected.x1[i] <=0 || Expected.x2[i] <=0){
        Expected.x3[i]  <- min(y1[i], y2[i])/2
        Expected.x1[i]  <- y1[i]-Expected.x3[i]
        Expected.x2[i]  <- y2[i]-Expected.x3[i]
        warning("negative expected x1 or x2 at ", i, " obs at ", j, " iteration", "\n")
      }
      den[i]              <- expected.latent.res$denominator
    }

    loglike.new    <- sum(den)
    loglike[j]     <- loglike.new
    loglike.diff <- abs(loglike.new-loglike.current)/(1+abs(loglike.new))
    if (verbose){
      cat(c("iter:", j, "; current loglike:", round(loglike.new,4), "; loglike change:", round(loglike.diff,6), "\n"))
    }

    #--------
    # M step
    #--------
    Q.b.function   <- function(coef){
      q.b.res <- sum((alpha1.current + alpha2.current + alpha3.current) * log(exp(coef %*% t(Model.Matrix.4)))) -
        sum(exp(coef %*% t(Model.Matrix.4)) * (Expected.x1 + Expected.x2 + Expected.x3))
      return(q.b.res)
    }
    beta.temp      <- (alpha1.current+alpha2.current+alpha3.current) / (Expected.x1+Expected.x2+Expected.x3)
    m4             <- stats::glm(beta.temp~Model.Matrix.4[,-1], family=Gamma(link="log"))
    coef4.new      <- stats::optim(par     = as.vector(m4$coefficients),
                            fn      = Q.b.function,
                            gr      = NULL,
                            hessian = TRUE,
                            control = list(fnscale=-1))$par
    beta.new       <- as.vector(exp(coef4.new%*%t(Model.Matrix.4)))

    alpha1.rootfun <- function(alpha1.var){
      sum(log(beta.new)) - digamma(alpha1.var)*n + sum(Expected.logx1)
    }
    alpha1.new <- stats::uniroot(alpha1.rootfun,
                          lower=sqrt(.Machine$double.eps),
                          upper=100000,
                          tol = sqrt(.Machine$double.xmin),
                          check.conv = TRUE)$root

    alpha2.rootfun <- function(alpha2.var){
      sum(log(beta.new)) - digamma(alpha2.var)*n + sum(Expected.logx2)
    }
    alpha2.new <- stats::uniroot(alpha2.rootfun,
                          lower=sqrt(.Machine$double.eps),
                          upper=100000,
                          tol = sqrt(.Machine$double.xmin),
                          check.conv = TRUE)$root
    alpha3.rootfun <- function(alpha3.var){
      sum(log(beta.new)) - digamma(alpha3.var)*n + sum(Expected.logx3)
    }
    alpha3.new <- stats::uniroot(alpha3.rootfun,
                          lower=sqrt(.Machine$double.eps),
                          upper=100000,
                          tol = sqrt(.Machine$double.xmin),
                          check.conv = TRUE)$root
    if (verbose){
      cat(c("alpha1:", round(alpha1.new,4), "\n"))
      cat(c("alpha2:", round(alpha2.new,4), "\n"))
      cat(c("alpha3:", round(alpha3.new,4), "\n"))
      cat(c("beta:",   round(summary(beta.new),4),  "\n"))
    }

    coef4.current  <- coef4.new
    loglike.current<- loglike.new
    alpha1.current <- alpha1.new
    alpha2.current <- alpha2.new
    alpha3.current <- alpha3.new
    beta.current   <- beta.new
    j              <- j+1
  }

  if (!all(loglike[1:(j-1)] == cummax(loglike[1:(j-1)]))){
    cat("Loglikelihood is not monotonically increasing!", "\n")
  }

  noparams<- m4$rank + 3
  AIC     <- -2*loglike[j-1] + noparams * 2
  BIC     <- -2*loglike[j-1] + noparams * log(n)

  # fitted values
  y1.fitted <- (alpha1.current + alpha3.current) / beta.current
  y2.fitted <- (alpha2.current + alpha3.current) / beta.current
  y1.residual <- y1 - y1.fitted
  y2.residual <- y2 - y2.fitted

  # formula
  f4.formula <- formula(paste("", as.character(f4.new)[3], sep="~"))

  result<-list(modelName   = "IE",
               coefficients= coef4.new,
               alpha1      = alpha1.current,
               alpha2      = alpha2.current,
               alpha3      = alpha3.current,
               beta        = beta.current,
               fitted.values= data.frame(y1=y1.fitted,
                                         y2=y2.fitted),
               residuals   = cbind(y1.residual, y2.residual),
               loglike     = loglike[j-1],
               ll          = loglike[1:(j-1)],
               df          = noparams,
               AIC         = AIC,
               BIC         = BIC,
               y           = cbind(y1, y2),
               n           = n,
               Model.Matrix= Model.Matrix.4,
               formula     = f4.formula,
               call        = tempcall,
               iter        = j-1)
  options(warn=0)
  structure(result, class = c('BGR'))
}

#' @export
print.BGR <- function (x, ...){
  txt <- paste0("'", class(x)[1], "' model of type '", x$modelName, "'")
  cat(txt, "\n")
  cat("\n")
  cat("Available components:\n")
  print(names(x))
  invisible(x)
}



