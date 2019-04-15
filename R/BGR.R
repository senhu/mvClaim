#' Bivariate gamma regression
#'
#' Bivariate gamma regression models for three model types: EE, EI, IE. Models are estimated by EM algorithm.
#'
#' @param modelName A character string indicating which model to be fitted. Need to be one of "EE", "EI", "IE".
#' @param response A vector of character strings indicating which variables in the data are treated as response or dependent variables.
#' @param data A matrix or data frame of observations. Categorical variables are allowed as covariates.
#' @param l1 A regression formula for the \eqn{latex}{\alpha_1} parameter in the bivariate gamma distribution. Note that, depending on the model type, might not be necessary to provide it.
#' @param l2 A regression formula for the \eqn{latex}{\alpha_2} parameter in the bivariate gamma distribution. Note that, depending on the model type, might not be necessary to provide it.
#' @param l3 A regression formula for the \eqn{latex}{\alpha_3} parameter in the bivariate gamma distribution. Note that, depending on the model type, might not be necessary to provide it.
#' @param l4 A regression formula for the \eqn{latex}{\beta} parameter in the bivariate gamma distribution. Note that, depending on the model type, might not be necessary to provide it.
#' @param expo Exposure variable
#' @param maxit A parameter that controls the number of maximum iteration in the EM algorithm. The default is 100.
#' @param tol A parameter that controls the convergence tolerance in the EM algorithm. The default is 1e-5.
#' @param Aitken A parameter that contols whether Aitken acceleration is used.
#' @param verbose A logical controlling whether estimations in each EM iteration are shown in the fitting process. The default is TRUE.
#'
#' @return An object of class \code{BGR} providing the estimation results.
#'   The details of the output components are:
#'   \item{call}{The matched call.}
#'   \item{coef}{The estimated coefficients.}
#'   \item{alpha1}{The estimated alpha1 values.}
#'   \item{alpha2}{The estimated alpha2 values.}
#'   \item{alpha3}{The estimated alpha3 values.}
#'   \item{beta}{The estimated beta values.}
#'   \item{trace}{All estimated coefficients and alpha, beta values in the EM algorithm.}
#'   \item{fitted}{The fitted values of the regression.}
#'   \item{loglike}{The final estimated maximum log-likelihood value.}
#'   \item{llseq}{The sequence of log-likelihood values in the EM algorithm fitting process.}
#'   \item{param.number}{Number of estimated parameters.}
#'   \item{AIC}{AIC values.}
#'   \item{BIC}{BIC values.}
#'   \item{iter}{Total iteration numbers.}
#'   \item{formula}{The formulas used in the regression.}
#'   \item{y}{The input response data.}
#'   \item{Model.Matrix}{The used model matrix for each regression formula.}
#'
#' @examples
#' m1 <- BGR(modelName = "EE",
#'           data   = regsim,
#'           response=c("y1","y2"),
#'           l1     = ~ w1 + w2,
#'           l2     = ~ w2 + w3,
#'           l3     = ~ w1 + w2 + w3,
#'           l4     = ~ w1 + w2 + w3,
#'           expo   = NULL,
#'           maxit  = 1000,
#'           tol    = 1e-6,
#'           Aitken = FALSE,
#'           verbose= TRUE)


BGR <- function(modelName = c("EE","EI","IE"),
                response,
                data,
                l1,
                l2,
                l3,
                l4,
                expo=NULL,
                maxit=100,
                tol=1e-5,
                Aitken = TRUE,
                verbose=TRUE){
  if (modelName == "EE"){
    # BGR.EE: all alpha's and beta are regressed on covariates
    if (is.null(l1)){ stop("l1 must be supplied for EE model type") }
    if (is.null(l2)){ stop("l2 must be supplied for EE model type") }
    if (is.null(l3)){ stop("l3 must be supplied for EE model type") }
    if (is.null(l4)){ stop("l4 must be supplied for EE model type") }
    return(BGR.EE(data=data, response=response,
                  l1=l1, l2=l2, l3=l3, l4=l4,
                  expo=expo, maxit=maxit, tol=tol,
                  Aitken = Aitken, verbose=verbose))
  }
  if (modelName == "EI"){
    # BGR.EI: same beta for all observations
    if (is.null(l1)){ stop("l1 must be supplied for EI model type") }
    if (is.null(l2)){ stop("l2 must be supplied for EI model type") }
    if (is.null(l3)){ stop("l3 must be supplied for EI model type") }
    return(BGR.EI(data=data, response=response,
                  l1=l1, l2=l2, l3=l3,
                  expo=expo, maxit=maxit, tol=tol,
                  Aitken=Aitken, verbose=verbose))
  }
  if (modelName == "IE"){
    # BGR.IE: only beta is regressed on its covariates, not alpha
    if (is.null(l4)){ stop("l4 must be supplied for IE model type") }
    return(BGR.IE(data=data, response=response,
                  l4=l4,
                  expo=expo, maxit=maxit, tol=tol,
                  Aitken=Aitken, verbose=verbose))
  }
}


#' @rdname BGR

BGR.EE <- function(data,
                   response,
                   l1,
                   l2,
                   l3,
                   l4,
                   expo   = NULL,
                   maxit  = 100,
                   tol    = 1e-5,
                   Aitken = FALSE,
                   verbose= TRUE){
  #------------------------------
  options(warn=-1)
  tempcall <- as.call( c(expression(BGR), list(data    = substitute(data),
                                               response= response,
                                               l1      = l1,
                                               l2      = l2,
                                               l3      = l3,
                                               l4      = l4,
                                               expo    = expo,
                                               maxit   = maxit,
                                               tol     = tol,
                                               Aitken  = Aitken,
                                               verbose = verbose) ) )
  #------------------------------
  namey1 <- response[1]
  namey2 <- response[2]
  y1     <- data[,names(data)==namey1]
  y2     <- data[,names(data)==namey2]
  n      <- length(y1)

  if (!is.null(expo)) {weight.var<-data[,names(data)==expo]} else {weight.var<-rep(1, n)}
  data$weight.var <- weight.var # there should be better ways of doing this

  data1 <- data;
  data1 <- data1[, names(data1)!=namey1];
  data1 <- data1[, names(data1)!=namey2]
  datap <- data1

  if (!is.null(l1)){
    if (as.character(l1[2])=="."){
      l1.new <- formula( paste( "x1", paste( names(data1),'',collapse='+',sep='' ), sep='~'))
    } else {l1.new <- formula(paste("x1", as.character(l1)[2], sep="~"))}
  } else stop("l1 cannot be NULL, must be supplied")
  if (!is.null(l2)){
    if (as.character(l2[2])=="."){l2.new <- formula( paste( "x2", paste( names(data1),'',collapse='+',sep='' ), sep='~'))}
    else {l2.new <- formula(paste("x2", as.character(l2)[2], sep="~"))}
  } else stop("l2 cannot be NULL, must be supplied")
  if (!is.null(l3)){
    if (as.character(l3[2])=="."){l3.new <- formula( paste( "x3", paste( names(data1), "", collapse = "+", sep = ""), sep="~"))}
    else {l3.new <- formula(paste("x3", as.character(l3)[2], sep="~"))}
  } else stop("l3 cannot be NULL, must be supplied")
  if (!is.null(l4)){
    if (as.character(l4[2])=="."){l4.new <- formula( paste( "be", paste( names(data1), "", collapse = "+", sep = ""), sep="~"))}
    else {l4.new <- formula(paste("be", as.character(l4)[2], sep="~"))}
  } else stop("l4 cannot be NULL, must be supplied")

  #-----------------------------
  # starting values
  x3.start <- rep(0, n)
  #set.seed(10)
  for (i in seq_len(n)){
    x3.start[i] <- runif(1, min=0, max=min(y1[i], y2[i]))
  }
  x1.start         <- y1 - x3.start
  x2.start         <- y2 - x3.start
  data.start       <- data1
  data.start$x1    <- x1.start
  data.start$x2    <- x2.start
  data.start$x3    <- x3.start
  start.temp.glm.1 <- glm(l1.new, data=data.start, family=Gamma(link="log"))
  start.temp.glm.2 <- glm(l2.new, data=data.start, family=Gamma(link="log"))
  start.temp.glm.3 <- glm(l3.new, data=data.start, family=Gamma(link="log"))

  Model.Matrix.1   <- model.matrix(start.temp.glm.1)
  Model.Matrix.2   <- model.matrix(start.temp.glm.2)
  Model.Matrix.3   <- model.matrix(start.temp.glm.3)
  coef1.current    <- start.temp.glm.1$coefficient
  coef2.current    <- start.temp.glm.2$coefficient
  coef3.current    <- start.temp.glm.3$coefficient

  alpha1.current   <- rep(1/summary(start.temp.glm.1)$dispersion, n) #rep(MASS::gamma.shape(start.temp.glm.1)$alpha, n)
  alpha2.current   <- rep(1/summary(start.temp.glm.2)$dispersion, n) #rep(MASS::gamma.shape(start.temp.glm.2)$alpha, n)
  alpha3.current   <- rep(1/summary(start.temp.glm.3)$dispersion, n) #rep(MASS::gamma.shape(start.temp.glm.3)$alpha, n)

  beta.current.1   <- 1/summary(start.temp.glm.1)$dispersion / predict(start.temp.glm.1, type="response")#, dispersion = 1/MASS::gamma.shape(start.temp.glm.1) )
  beta.current.2   <- 1/summary(start.temp.glm.2)$dispersion / predict(start.temp.glm.2, type="response")#, dispersion = 1/MASS::gamma.shape(start.temp.glm.2) )
  beta.current.3   <- 1/summary(start.temp.glm.3)$dispersion / predict(start.temp.glm.3, type="response")#, dispersion = 1/MASS::gamma.shape(start.temp.glm.3) )
  beta.current     <- rowMeans(cbind(beta.current.1, beta.current.2, beta.current.3))

  data.start$be    <- beta.current
  start.temp.glm.4 <- glm(l4.new, data=data.start, family=Gamma(link="log"))
  coef4.current    <- start.temp.glm.4$coefficient
  Model.Matrix.4   <- model.matrix(start.temp.glm.4)

  # p1 <- dim(Model.Matrix.1)[2]
  # p2 <- dim(Model.Matrix.2)[2]
  # p3 <- dim(Model.Matrix.3)[2]
  # p4 <- dim(Model.Matrix.4)[2]
  #--------------------------
  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)
  loglike.store   <- rep(-.Machine$double.xmax, 3)
  #if (is.null(l3)) stop("l3 is NULL, Double Gamma Regression")

  coef1.trace  <- coef1.current
  coef2.trace  <- coef2.current
  coef3.trace  <- coef3.current
  coef4.trace  <- coef4.current
  alpha1.trace <- alpha1.current
  alpha2.trace <- alpha2.current
  alpha3.trace <- alpha3.current
  beta.trace   <- beta.current
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
      expected.latent.res <- expected.latent(c(y1[i],y2[i]),
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
    loglike.store  <- c(loglike.store[-1], loglike.new)
    if (Aitken){
      loglike.diff <- abs( MoEClust::MoE_aitken(loglike.store)$linf - loglike.current )
    } else {
      loglike.diff <- abs(loglike.new-loglike.current)/(1+abs(loglike.new))
    }
    if (verbose){
      cat(c("iter:", j, "; current loglike:", round(loglike.new,4), "; loglike change:", round(loglike.diff,6), "\n"))
    }

    #--------
    # M step
    #--------

    beta.temp      <- (alpha1.current + alpha2.current + alpha3.current) / (Expected.x1 + Expected.x2 + Expected.x3)
    m4             <- glm(beta.temp ~ Model.Matrix.4[,-1], family=Gamma(link="log"))
    Q.b.function   <- function(coef){
      q.b.res <- sum((alpha1.current + alpha2.current + alpha3.current) * log(exp(coef %*% t(Model.Matrix.4)))) -
        sum(exp(coef %*% t(Model.Matrix.4)) * (Expected.x1 + Expected.x2 + Expected.x3))
      return(q.b.res)
    }
    coef4.optim    <- optim(par     = as.vector(m4$coefficients),
                            fn      = Q.b.function,
                            gr      = NULL,
                            hessian = TRUE,
                            control = list(fnscale=-1, maxit=1e+5))
    if (coef4.optim$convergence!=0) cat("optim coef. beta not converged")
    coef4.new      <- coef4.optim$par
    beta.new       <- as.vector(exp(coef4.new%*%t(Model.Matrix.4)))

    m1             <- glm(alpha1.current~Model.Matrix.1[,-1], family=Gamma(link="log"))
    Q.function.a1  <- function(coef){
      q.res <- sum( exp(coef %*% t(Model.Matrix.1)) * log(beta.new) ) -
        sum( lgamma(exp(coef%*%t(Model.Matrix.1))) ) +
        sum( exp(coef %*% t(Model.Matrix.1)) * Expected.logx1 )
      return(q.res)
    }
    coef1.optim    <- optim(par     = as.vector(m1$coefficients),
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
    m2             <- glm(alpha2.current~Model.Matrix.2[,-1], family=Gamma(link="log"))
    coef2.optim    <- optim(par     = as.vector(m2$coefficients),
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
    m3             <- glm(alpha3.current~Model.Matrix.3[,-1], family=Gamma(link="log"))
    coef3.optim    <- optim(par     = as.vector(m3$coefficients),
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
    alpha1.trace   <- cbind(alpha1.trace, alpha1.new)
    alpha2.trace   <- cbind(alpha2.trace, alpha2.new)
    alpha3.trace   <- cbind(alpha3.trace, alpha3.new)
    beta.trace     <- cbind(beta.trace, beta.new)
    coef1.trace    <- cbind(coef1.trace, coef1.new)
    coef2.trace    <- cbind(coef2.trace, coef2.new)
    coef3.trace    <- cbind(coef3.trace, coef3.new)
    coef4.trace    <- cbind(coef4.trace, coef4.new)
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
  l1.formula <- formula(paste("", as.character(l1.new)[3], sep="~"))
  l2.formula <- formula(paste("", as.character(l2.new)[3], sep="~"))
  l3.formula <- formula(paste("", as.character(l3.new)[3], sep="~"))
  l4.formula <- formula(paste("", as.character(l4.new)[3], sep="~"))

  result<-list(coef        =list(coef1.new, coef2.new, coef3.new),
               alpha1      =alpha1.current,
               alpha2      =alpha2.current,
               alpha3      =alpha3.current,
               beta        =beta.current,
               trace       =list(alpha1.trace=alpha1.trace,
                                 alpha2.trace=alpha2.trace,
                                 alpha3.trace=alpha3.trace,
                                 beta.trace  =beta.trace,
                                 coef1.trace  =coef1.trace,
                                 coef2.trace  =coef2.trace,
                                 coef3.trace  =coef3.trace,
                                 coef4.trace  =coef4.trace),
               fitted      =data.frame(y1=y1.fitted,
                                       y2=y2.fitted),
               loglike     =loglike[j-1],
               llseq       = loglike[1:(j-1)],

               param.number= noparams,
               AIC         = AIC,
               BIC         = BIC,

               y           = cbind(y1, y2),
               Model.Matrix= list(Model.Matrix.1, Model.Matrix.2, Model.Matrix.3),
               formula     = list(l1.formula, l2.formula, l3.formula, l4.formula),

               call        = tempcall,
               iter        = j-1)
  options(warn=0)
  class(result)<-c('BGR.EE', 'glm')
  return(result)
}

#' @rdname BGR

BGR.EI <- function(data,
                   response,
                   l1,
                   l2,
                   l3,
                   expo   = NULL,
                   maxit  = 100,
                   tol    = 1e-5,
                   Aitken = FALSE,
                   verbose= TRUE){
  #------------------------------
  options(warn=-1)
  tempcall <- as.call( c(expression(BGR), list(data    = substitute(data),
                                               response= response,
                                               l1      = l1,
                                               l2      = l2,
                                               l3      = l3,
                                               expo    = expo,
                                               maxit   = maxit,
                                               tol     = tol,
                                               Aitken  = Aitken,
                                               verbose = verbose) ) )
  #------------------------------
  namey1 <- response[1]
  namey2 <- response[2]
  y1     <- data[,names(data)==namey1]
  y2     <- data[,names(data)==namey2]
  n      <- length(y1)

  if (!is.null(expo)) {weight.var<-data[,names(data)==expo]} else {weight.var<-rep(1, n)}
  data$weight.var <- weight.var # there should be better ways of doing this

  if (!is.null(l1)){
    if (as.character(l1[2])=="."){
      l1.new <- formula( paste( "x1", paste( names(data1),'',collapse='+',sep='' ), sep='~'))
    } else {l1.new <- formula(paste("x1", as.character(l1)[2], sep="~"))}
  } else stop("l1 cannot be NULL, must be supplied")
  if (!is.null(l2)){
    if (as.character(l2[2])=="."){l2.new <- formula( paste( "x2", paste( names(data1),'',collapse='+',sep='' ), sep='~'))}
    else {l2.new <- formula(paste("x2", as.character(l2)[2], sep="~"))}
  } else stop("l2 cannot be NULL, must be supplied")
  if (!is.null(l3)){
    if (as.character(l3[2])=="."){l3.new <- formula( paste( "x3", paste( names(data1), "", collapse = "+", sep = ""), sep="~"))}
    else {l3.new <- formula(paste("x3", as.character(l3)[2], sep="~"))}
  } else stop("l3 cannot be NULL, must be supplied")

  data1 <- data;
  data1<-data1[, names(data1)!=namey1];
  data1<-data1[, names(data1)!=namey2]
  datap <- data1

  #-----------------------------
  # starting values
  x3.start <- rep(0, n)
  #set.seed(10)
  for (i in seq_len(n)){
    x3.start[i] <- runif(1, min=0, max=min(y1[i], y2[i]))
  }
  x1.start       <- y1 - x3.start
  x2.start       <- y2 - x3.start
  data.start     <- data1
  data.start$x1  <- x1.start
  data.start$x2  <- x2.start
  data.start$x3  <- x3.start
  start.temp.glm.1 <- glm(l1.new, data=data.start, family=Gamma(link="log"))
  start.temp.glm.2 <- glm(l2.new, data=data.start, family=Gamma(link="log"))
  start.temp.glm.3 <- glm(l3.new, data=data.start, family=Gamma(link="log"))
  Model.Matrix.1 <- stats::model.matrix(start.temp.glm.1)
  Model.Matrix.2 <- stats::model.matrix(start.temp.glm.2)
  Model.Matrix.3 <- stats::model.matrix(start.temp.glm.3)
  coef1.current  <- start.temp.glm.1$coefficient
  coef2.current  <- start.temp.glm.2$coefficient
  coef3.current  <- start.temp.glm.3$coefficient

  alpha1.current <- rep(1/summary(start.temp.glm.1)$dispersion, n) #rep(MASS::gamma.shape(start.temp.glm.1)$alpha, n)
  alpha2.current <- rep(1/summary(start.temp.glm.2)$dispersion, n) #rep(MASS::gamma.shape(start.temp.glm.2)$alpha, n)
  alpha3.current <- rep(1/summary(start.temp.glm.3)$dispersion, n) #rep(MASS::gamma.shape(start.temp.glm.3)$alpha, n)
  beta.current.1 <- 1/summary(start.temp.glm.1)$dispersion / predict(start.temp.glm.1, type="response")#, dispersion = 1/MASS::gamma.shape(start.temp.glm.1) )
  beta.current.2 <- 1/summary(start.temp.glm.2)$dispersion / predict(start.temp.glm.2, type="response")#, dispersion = 1/MASS::gamma.shape(start.temp.glm.2) )
  beta.current.3 <- 1/summary(start.temp.glm.3)$dispersion / predict(start.temp.glm.3, type="response")#, dispersion = 1/MASS::gamma.shape(start.temp.glm.3) )
  beta.current   <- mean(cbind(beta.current.1, beta.current.2, beta.current.3))
  #--------------------------

  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)
  loglike.store   <- rep(-.Machine$double.xmax, 3)
  #if (is.null(l3)) stop("l3 is NULL, Double Gamma Regression")

  coef1.trace  <- coef1.current
  coef2.trace  <- coef2.current
  coef3.trace  <- coef3.current
  alpha1.trace <- alpha1.current
  alpha2.trace <- alpha2.current
  alpha3.trace <- alpha3.current
  beta.trace   <- beta.current
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
      expected.latent.res <- expected.latent(c(y1[i],y2[i]),
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
    loglike.store  <- c(loglike.store[-1], loglike.new)
    if (Aitken){
      loglike.diff <- abs( MoEClust::MoE_aitken(loglike.store)$linf - loglike.current )
    } else {
      loglike.diff <- abs(loglike.new-loglike.current)/(1+abs(loglike.new))
    }
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

    m1             <- glm(alpha1.current~Model.Matrix.1[,-1], family=Gamma(link="log"))
    coef1.new      <- optim(par          = as.vector(m1$coefficients),
                            fn           = Q.function,
                            Model.Matrix = Model.Matrix.1,
                            Expected.logs= Expected.logx1,
                            gr           = NULL,
                            hessian      = TRUE,
                            control      = list(fnscale=-1))$par
    alpha1.new     <- as.vector(exp(coef1.new%*%t(Model.Matrix.1)))

    m2             <- glm(alpha2.current~Model.Matrix.2[,-1], family=Gamma(link="log"))
    coef2.new      <- optim(par          = as.vector(m2$coefficients),
                            fn           = Q.function,
                            Model.Matrix = Model.Matrix.2,
                            Expected.logs= Expected.logx2,
                            gr           = NULL,
                            hessian      = TRUE,
                            control      = list(fnscale=-1))$par
    alpha2.new     <- as.vector(exp(coef2.new%*%t(Model.Matrix.2)))

    m3             <- glm(alpha3.current~Model.Matrix.3[,-1], family=Gamma(link="log"))
    coef3.new      <- optim(par          = as.vector(m3$coefficients),
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

    coef1.current  <- coef1.new
    coef2.current  <- coef2.new
    coef3.current  <- coef3.new

    loglike.current<- loglike.new

    alpha1.current <- alpha1.new
    alpha2.current <- alpha2.new
    alpha3.current <- alpha3.new
    beta.current   <- beta.new

    alpha1.trace   <- cbind(alpha1.trace, alpha1.new)
    alpha2.trace   <- cbind(alpha2.trace, alpha2.new)
    alpha3.trace   <- cbind(alpha3.trace, alpha3.new)
    beta.trace     <- c(beta.trace, beta.new)
    coef1.trace    <- cbind(coef1.trace, coef1.new)
    coef2.trace    <- cbind(coef2.trace, coef2.new)
    coef3.trace    <- cbind(coef3.trace, coef3.new)

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
  l1.formula <- formula(paste("", as.character(l1.new)[3], sep="~"))
  l2.formula <- formula(paste("", as.character(l2.new)[3], sep="~"))
  l3.formula <- formula(paste("", as.character(l3.new)[3], sep="~"))

  result<-list(coef        =list(coef1.new, coef2.new, coef3.new),
               alpha1      =alpha1.current,
               alpha2      =alpha2.current,
               alpha3      =alpha3.current,
               beta        =beta.current,
               trace       =list(alpha1.trace=alpha1.trace,
                                 alpha2.trace=alpha2.trace,
                                 alpha3.trace=alpha3.trace,
                                 beta.trace  =beta.trace,
                                 coef1.trace =coef1.trace,
                                 coef2.trace =coef2.trace,
                                 coef3.trace =coef3.trace),
               fitted      =data.frame(y1=y1.fitted,
                                       y2=y2.fitted),
               loglike     =loglike[j-1],
               llseq       = loglike[1:(j-1)],

               param.number=noparams,
               AIC         = AIC,
               BIC         = BIC,

               y           = cbind(y1, y2),
               Model.Matrix= list(Model.Matrix.1, Model.Matrix.2, Model.Matrix.3),
               formula     = list(l1.formula, l2.formula, l3.formula),

               call        = tempcall,
               iter        = j-1)
  options(warn=0)
  class(result)<-c('BGR.EI', 'glm')
  return(result)
}

#' @rdname BGR

BGR.IE <- function(data,
                   response,
                   l4,
                   expo   = NULL,
                   maxit  = 100,
                   tol    = 1e-5,
                   Aitken = FALSE,
                   verbose= TRUE){
  #------------------------------
  options(warn=-1)
  tempcall <- as.call( c(expression(BGR), list(data    = substitute(data),
                                               response= response,
                                               l4      = l4,
                                               expo    = expo,
                                               maxit   = maxit,
                                               tol     = tol,
                                               Aitken  = Aitken,
                                               verbose = verbose) ) )
  #------------------------------
  namey1 <- response[1]
  namey2 <- response[2]
  y1     <- data[,names(data)==namey1]
  y2     <- data[,names(data)==namey2]
  n      <- length(y1)

  if (!is.null(expo)) {weight.var<-data[,names(data)==expo]} else {weight.var<-rep(1, n)}
  data$weight.var <- weight.var

  data1 <- data;
  data1 <- data1[,names(data1)!=namey1];
  data1 <- data1[,names(data1)!=namey2]
  datap <- data1

  if (!is.null(l4)){
    if (as.character(l4[2])=="."){l4.new <- formula( paste( "be", paste( names(data1), "", collapse = "+", sep = ""), sep="~"))}
    else {l4.new <- formula(paste("be", as.character(l4)[2], sep="~"))}
  } else stop("l4 cannot be NULL, must be supplied")

  #-----------------------------
  # starting values
  x3.start <- rep(0, n) # set.seed(10)
  for (i in seq_len(n)){
    x3.start[i] <- runif(1, min=0, max=min(y1[i], y2[i]))
  }
  x1.start         <- y1 - x3.start
  x2.start         <- y2 - x3.start

  alpha1.current   <- MASS::fitdistr(x1.start, "gamma")$estimate[1]
  alpha2.current   <- MASS::fitdistr(x1.start, "gamma")$estimate[1]
  alpha3.current   <- MASS::fitdistr(x1.start, "gamma")$estimate[1]
  beta.current     <- rowMeans(cbind(alpha1.current/x1.start, alpha2.current/x2.start, alpha3.current/x3.start))

  data.start       <- data1
  data.start$be    <- beta.current
  start.temp.glm   <- glm(l4.new, data=data.start, family=Gamma(link="log"))
  coef4.current    <- start.temp.glm$coefficient
  Model.Matrix.4   <- model.matrix(start.temp.glm)

  #--------------------------

  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)
  loglike.store   <- rep(-.Machine$double.xmax, 3)

  coef4.trace  <- coef4.current
  alpha1.trace <- alpha1.current
  alpha2.trace <- alpha2.current
  alpha3.trace <- alpha3.current
  beta.trace   <- beta.current
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
      expected.latent.res <- expected.latent(c(y1[i],y2[i]),
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
    loglike.store  <- c(loglike.store[-1], loglike.new)
    if (Aitken){
      loglike.diff <- abs( MoEClust::MoE_aitken(loglike.store)$linf - loglike.current )
    } else {
      loglike.diff <- abs(loglike.new-loglike.current)/(1+abs(loglike.new))
    }
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
    m4             <- glm(beta.temp~Model.Matrix.4[,-1], family=Gamma(link="log"))
    coef4.new      <- optim(par     = as.vector(m4$coefficients),
                            fn      = Q.b.function,
                            gr      = NULL,
                            hessian = TRUE,
                            control = list(fnscale=-1))$par
    beta.new       <- as.vector(exp(coef4.new%*%t(Model.Matrix.4)))

    alpha1.rootfun <- function(alpha1.var){
      sum(log(beta.new)) - digamma(alpha1.var)*n + sum(Expected.logx1)
    }
    alpha1.new <- uniroot(alpha1.rootfun,
                          lower=sqrt(.Machine$double.eps),
                          upper=100000,
                          tol = sqrt(.Machine$double.xmin),
                          check.conv = TRUE)$root

    alpha2.rootfun <- function(alpha2.var){
      sum(log(beta.new)) - digamma(alpha2.var)*n + sum(Expected.logx2)
    }
    alpha2.new <- uniroot(alpha2.rootfun,
                          lower=sqrt(.Machine$double.eps),
                          upper=100000,
                          tol = sqrt(.Machine$double.xmin),
                          check.conv = TRUE)$root
    alpha3.rootfun <- function(alpha3.var){
      sum(log(beta.new)) - digamma(alpha3.var)*n + sum(Expected.logx3)
    }
    alpha3.new <- uniroot(alpha3.rootfun,
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
    alpha1.trace   <- cbind(alpha1.trace, alpha1.new)
    alpha2.trace   <- cbind(alpha2.trace, alpha2.new)
    alpha3.trace   <- cbind(alpha3.trace, alpha3.new)
    beta.trace     <- cbind(beta.trace, beta.new)
    coef4.trace    <- cbind(coef4.trace, coef4.new)
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
  l4.formula <- formula(paste("", as.character(l4.new)[3], sep="~"))

  result<-list(coef        =coef4.new,
               alpha1      =alpha1.current,
               alpha2      =alpha2.current,
               alpha3      =alpha3.current,
               beta        =beta.current,
               trace       =list(alpha1.trace=alpha1.trace,
                                 alpha2.trace=alpha2.trace,
                                 alpha3.trace=alpha3.trace,
                                 beta.trace  =beta.trace,
                                 coef.trace  =coef4.trace),
               fitted      =data.frame(y1=y1.fitted,
                                       y2=y2.fitted),
               loglike     =loglike[j-1],
               llseq       = loglike[1:(j-1)],

               param.number=noparams,
               AIC         = AIC,
               BIC         = BIC,

               y           = cbind(y1, y2),
               Model.Matrix= list(Model.Matrix.4),
               formula     = list(l4.formula),

               call        = tempcall,
               iter        = j-1)
  options(warn=0)
  class(result)<-c('BGR.IE', 'glm')
  return(result)
}





