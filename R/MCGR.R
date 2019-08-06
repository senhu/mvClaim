#' Mixture of Copula Regressions with Gamma Marginal Distributions
#'
#' Mixture of copula regressions, or model-based clustering with copulas and covariates. Models are estimated by EM algorithm.
#'
#' @param copula A specified copula, from the \code{copula} package. Must be a list of copulas or a list of one copula, matching G
#' @param f1 A regression formula for the first marginal GLM.
#' @param f2 A regression formula for the second marginal GLM.
#' @param data A matrix or data frame of observations. Categorical variables are allowed as covariates.
#' @param G An interger vector specifying the number of mixture components (clusters).
#' @param gating Specifies the gating network, can be "C", "E" or a regression formula.
#' @param initialization Specifies initialization method for EM algorithm. The default is "\code{mclust}".
#' @param maxit A parameter that controls the number of maximum iteration in the EM algorithm. The default is 100.
#' @param tol A parameter that controls the convergence tolerance in the EM algorithm. The default is 1e-5.
#' @param verbose A logical controlling whether estimations in each EM iteration are shown in the fitting process. The default is TRUE.
#'
#' @return An object of class \code{BGR} providing the estimation results.
#'   The details of the output components are:
#'   \item{call}{The matched call.}
#'   \item{coefficients}{The estimated coefficients in the expert and gating networks, if exists.}
#'   \item{copula.param}{The estimated copula parameters.}
#'   \item{copula}{The specified copula.}
#'   \item{theta}{All estimated parameter values, in a vector format.}
#'   \item{pro}{A vector whose g-th component is the mixing proportion for the g-th component of the mixture model.}
#'   \item{z}{A matrix whose [i,g]-th entry is the probability that observation i in the data belongs to the g-th group.}
#'   \item{class}{The classification corresponding to z.}
#'   \item{fitted.values}{The fitted values of the regression.}
#'   \item{residuals}{The residuals of the regression}
#'   \item{loglike}{The final estimated maximum log-likelihood value.}
#'   \item{ll}{The sequence of log-likelihood values in the EM algorithm fitting process.}
#'   \item{df}{The number of estimated parameters.}
#'   \item{AIC}{AIC values.}
#'   \item{BIC}{BIC values.}
#'   \item{iter}{Total iteration numbers.}
#'   \item{formula}{The formulas used in the regression.}
#'   \item{y}{The input response data.}
#'   \item{n}{The number of observations in the data.}
#'   \item{gating.model}{The binomial/multinomial regression model in the gating network.}
#'   \item{Model.Matrix}{The used model matrix for each regression formula.}
#'
#' @examples
#' \donttest{
#' data("simdat.mcgr")
#' res <- MCGR(copula = list(copula::gumbelCopula(dim=2),
#'                           copula::frankCopula(dim=2)),
#'             f1 = y1 ~ x1+x2,
#'             f2 = y2 ~ x1+x2,
#'             G = 2,
#'             gating = "C",
#'             data = simdat.mcgr,
#'             verbose = FALSE)
#' }
#'
#' @export MCGR

MCGR <- function(copula,
                 f1,
                 f2,
                 G,
                 data,
                 gating,
                 initialization = "mclust",
                 maxit   = 100,
                 tol     = 1e-5,
                 verbose = TRUE){
  tempcall <- as.call(c(expression(MCGR),list(copula = copula,
                                              f1     = f1,
                                              f2     = f2,
                                              data   = substitute(data),
                                              G      = G,
                                              gating = gating,
                                              initialization = initialization,
                                              maxit  = maxit,
                                              tol    = tol,
                                              verbose= verbose)))

  n <- nrow(data)
  namey1 <- as.character(f1[2])
  namey2 <- as.character(f2[2])
  y1     <- data[,names(data)==namey1]
  y2     <- data[,names(data)==namey2]
  y      <- cbind(y1,y2)

  newdata <- data
  newdata <- newdata[,names(newdata)!=namey1]
  newdata <- newdata[,names(newdata)!=namey2]

  if (!is.null(f2) && class(f1)=="formula" && as.character(f1[3])=="."){
    f1 <- stats::formula( paste( "y1", paste( names(newdata),'',collapse='+',sep='' ), sep='~'))
  }
  if (!is.null(f2) && class(f2)=="formula" && as.character(f2[3])=="."){
    f2 <- stats::formula( paste( "y2", paste( names(newdata),'',collapse='+',sep='' ), sep='~'))
  }
  if (gating!="C" && gating != "E"){
    if (as.character(gating[2])=="."){gating.formula <- stats::formula( paste( "pzy", paste( names(newdata), "", collapse = "+", sep = ""), sep="~"))}
    else {gating.formula <- stats::formula(paste("pzy", as.character(gating)[2], sep="~"))}
  }
  model.matrix.1 <- stats::model.matrix(f1[-2], data = newdata)
  model.matrix.2 <- stats::model.matrix(f2[-2], data = newdata)
  l1 <- ncol(model.matrix.1)+1
  l2 <- ncol(model.matrix.2)+1

  if (length(copula) == G){
    coplist <- copula
  } else if (G>1 && length(copula)==1){
    coplist <- rep(list(copula),G)
  } else {stop("Provided copulas do not match G value")}

  if (initialization == "mclust"){
    if (G>1){
      initial.mclust <- mclust::Mclust(G=G, data=y)
      z.init         <- initial.mclust$z
      if (gating == "C"){
        tau.current        <- initial.mclust$parameters$pro
      } else if (gating == "E"){
        tau.current      <- rep(1/G, G)
      } else {
        newdata$pzy      <- z.init
        mp             <- nnet::multinom(formula = gating.formula, data = newdata, trace=FALSE)  # reltol
        tau.current    <- stats::predict(mp, type="probs") }
      index <- initial.mclust$classification
      m1 <- stats::glm(f1,family=Gamma(link="log"),data=data)
      m2 <- stats::glm(f2,family=Gamma(link="log"),data=data)
      theta.current  <- NULL
      for (gg in seq_len(G)){
        a.temp <- tryCatch({
          copula::fitCopula(coplist[[gg]],copula::pobs(y),
                            method = "mpl", optim.method="Nelder-Mead")@estimate
        }, error = function(err){
          err.temp <- copula::fitCopula(coplist[[gg]],copula::pobs(y[index==gg,]),
                                        method = "itau")@estimate
          return(err.temp)
        })
        theta.current <- append(theta.current,
                                list(c(m1$coefficients,
                                       log(MASS::gamma.shape(m1)$alpha),
                                       m2$coefficients,
                                       log(MASS::gamma.shape(m2)$alpha),
                                       a.temp)))
      }
    }
    if (G==1){
      z.init         <- rep(1,n)
      if (gating == "C" || gating == "E"){
        tau.current      <- 1
      } else {
        tau.current    <- rep(1,n) }
      index <- rep(1,n)
      m1 <- stats::glm(f1,family=Gamma(link="log"),data=data)
      m2 <- stats::glm(f2,family=Gamma(link="log"),data=data)
      a.temp <- tryCatch({
        copula::fitCopula(coplist[[1]],copula::pobs(y),
                          method = "mpl", optim.method="Nelder-Mead")@estimate
      }, error = function(err){
        err.temp <- copula::fitCopula(coplist[[1]],copula::pobs(y),
                                      method = "itau")@estimate
        return(err.temp)
      })
      theta.current <- list(c(m1$coefficients,
                              log(MASS::gamma.shape(m1)$alpha),
                              m2$coefficients,
                              log(MASS::gamma.shape(m2)$alpha),
                              a.temp))
    }
  } #end of initialization
  if (class(initialization) == "Mclust"){
    if (G>1){
      z.init         <- initialization$z
      if (gating == "C"){
        tau.current        <- initialization$parameters$pro
      } else if (gating == "E"){
        tau.current      <- rep(1/G, G)
      } else {
        newdata$pzy      <- z.init
        mp             <- nnet::multinom(formula = gating.formula, data = newdata, trace=FALSE)  # reltol
        tau.current    <- stats::predict(mp, type="probs") }
      index <- initialization$classification
      m1 <- stats::glm(f1,family=Gamma(link="log"),data=data)
      m2 <- stats::glm(f2,family=Gamma(link="log"),data=data)
      theta.current  <- NULL
      for (gg in seq_len(G)){
        a.temp <- tryCatch({
          copula::fitCopula(coplist[[gg]],copula::pobs(y),
                            method = "mpl", optim.method="Nelder-Mead")@estimate
        }, error = function(err){
          err.temp <- copula::fitCopula(coplist[[gg]],copula::pobs(y[index==gg,]),
                                        method = "itau")@estimate
          return(err.temp)
        })
        theta.current <- append(theta.current,
                                list(c(m1$coefficients,
                                       log(MASS::gamma.shape(m1)$alpha),
                                       m2$coefficients,
                                       log(MASS::gamma.shape(m2)$alpha),
                                       a.temp)))
      }
    }
    if (G==1){
      z.init         <- rep(1,n)
      if (gating == "C" || gating == "E"){
        tau.current      <- 1
      } else {
        tau.current    <- rep(1,n) }
      index <- rep(1,n)
      m1 <- stats::glm(f1,family=Gamma(link="log"),data=data)
      m2 <- stats::glm(f2,family=Gamma(link="log"),data=data)
      a.temp <- tryCatch({
        copula::fitCopula(coplist[[1]],copula::pobs(y),
                          method = "mpl", optim.method="Nelder-Mead")@estimate
      }, error = function(err){
        err.temp <- copula::fitCopula(coplist[[1]], copula::pobs(y),
                                      method = "itau")@estimate
        return(err.temp)
      })
      theta.current <- list(c(m1$coefficients,
                              log(MASS::gamma.shape(m1)$alpha),
                              m2$coefficients,
                              log(MASS::gamma.shape(m2)$alpha),
                              a.temp))
    }
  } #end of initialization

  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)
  j               <- 1
  while ( (loglike.diff > tol) && (j <= maxit) ) {
    #-------------
    # E-step
    #-------------
    logden <- matrix(0,nrow=n, ncol=G)
    for (i in seq_len(n)){
      for (g in seq_len(G)){
        logden[i,g] <- thelogden.gamma(theta=theta.current[[g]],
                                       y=y[i,],
                                       xmat=list(model.matrix.1[i,],
                                            model.matrix.2[i,]),
                                       copula=coplist[[g]])
        #print(c(i,g))
      }
    }
    if (gating=="C" || gating=="E"){
      newlogden   <- sweep(logden, 2, log(tau.current), FUN="+") } else{
        newlogden <- logden + log(tau.current)  }
    logdenom   <- matrixStats::rowLogSumExps(newlogden)
    z.new   <- exp(sweep(newlogden, 1, logdenom, FUN="-"))

    loglike.new    <- sum(logdenom)
    loglike[j]     <- loglike.new
    loglike.diff <- abs(loglike.new-loglike.current)/(1+abs(loglike.new))
    if (verbose){
      cat(c("iter:", j, "; current loglike:", round(loglike.new,4), "; loglike change:", round(loglike.diff,4), "\n"))
    }
    #-------------
    # M-step
    #-------------
    if (gating == "C"){
      tau.new        <- colSums(z.new)/n
      coef.p.new     <- NULL
      mp             <- NULL } else if (gating == "E"){
        tau.new      <- rep(1/G, G)
        coef.p.new   <- NULL
        mp           <- NULL } else {
          newdata$pzy<- z.new
          mp         <- nnet::multinom(formula = gating.formula, data = newdata, trace=FALSE)  # reltol
          coef.p.new <- stats::coef(mp)
          tau.new    <- stats::predict(mp, type="probs") }

    # max all parameters together
    theta.new <- theta.current
    for (g in c(1:G)){
      fit <- stats::optim(theta.current[[g]], Q.function,
                          control = list(fnscale = -1, trace = 0),
                          y       = y,
                          xmat    = list(model.matrix.1, model.matrix.2),
                          copula  = coplist[[g]],
                          z.mat   = z.new[,g],
                          method  = "Nelder-Mead")
      #if (fit$convergence != 0) cat("optim not converged at G=",g, ". \n")
      theta.new[[g]] <- fit$par
    }

    if (verbose){
      if (is.matrix(tau.new)) { cat(c("tau:", round(colMeans(tau.new),4), "\n")) } else cat(c("tau:", round(tau.new,3), "\n"))
    }
    loglike.current<- loglike.new
    theta.current  <- theta.new
    tau.current    <- tau.new
    j              <- j+1
  } # end of EM

  if (!all(loglike[1:(j-1)] == cummax(loglike[1:(j-1)]))){
    cat("Loglikelihood is not strictly monotonically increasing!", "\n")
  }
  if (gating=="C"){
    noparams <- sum(sapply(theta.current, length))+(G-1)
  } else if (gating=="E"){
    noparams <- sum(sapply(theta.current, length))
  } else {
    noparams <- sapply(theta.current, length)+mp$rank
  }

  AIC<- -2*loglike[j-1] + noparams * 2
  BIC<- -2*loglike[j-1] + noparams * log(n)
  classification <- apply(z.new, 1, which.max)

  coef1 <- rep(list(NULL),G)
  coef2 <- rep(list(NULL),G)
  shape1 <- rep(0,G)
  shape2 <- rep(0,G)
  pcop <- rep(list(NULL),G)
  fitval1 <- matrix(0, nrow=n, ncol=G)
  fitval2 <- matrix(0, nrow=n, ncol=G)
  for (g in c(1:G)){
    finparam <- theta.current[[g]]
    coef1[[g]] <- finparam[1:(l1-1)]
    coef2[[g]] <- finparam[(l1+1):(l1+l2-1)]
    pcop[[g]] <- finparam[-(1:(l1 + l2))]
    fitval1[,g] <- exp(model.matrix.1%*%coef1[[g]])
    fitval2[,g] <- exp(model.matrix.2%*%coef2[[g]])
    shape1[g] <- exp(finparam[l1])
    shape2[g] <- exp(finparam[(l1+l2)])
    if (class(coplist[[g]])=="rotExplicitCopula"){
      coplist[[g]]@copula@parameters <- pcop[[g]]
    } else {coplist[[g]]@parameters <- pcop[[g]]}
    #coplist[[g]]@parameters <- pcop[[g]]
  }

  if (gating=="C" || gating=="E"){
    finfitval <- cbind(fitval1%*%tau.current, fitval2%*%tau.current)
  } else {
    finfitval <- cbind(rowSums(fitval1*tau.current), rowSums(fitval2*tau.current))
  }
  residuals <- (y - finfitval)

  f1.formula <- stats::formula(paste(namey1, as.character(f1)[3], sep="~"))
  f2.formula <- stats::formula(paste(namey2, as.character(f2)[3], sep="~"))
  if (gating == "C") {
    gating.formula <- "C"} else if (gating == "E") {
      gating.formula <- "E"} else {
        gating.formula <- formula(paste("", as.character(gating.formula)[3], sep="~")) }

  #  Calculation of output
  result<-list(coefficients = c(beta1=coef1,
                                shape1=list(shape1),
                                beta2=coef2,
                                shape2=list(shape2)),
               gating.coef  = coef.p.new,
               copula.param = pcop,
               copula       = coplist,
               theta        = theta.current,
               pro          = tau.current,
               z            = z.new,
               class        = classification,
               fitted.values= data.frame(finfitval),
               residuals    = residuals,
               loglike      = loglike.current,
               ll           = loglike[1:(j-1)],
               df           = noparams,
               AIC          = AIC,
               BIC          = BIC,
               n            = n,
               y            = cbind(y1, y2),
               Model.Matrix = list(model.matrix.1,
                                   model.matrix.2),
               formula      = list(f1 = f1.formula,
                                   f2 = f2.formula,
                                   gating = gating.formula),
               gating.model = mp,
               call         = tempcall,
               iter         = (j-1))
  options(warn=0)
  structure(result, class = 'MCGR')
}

loglik.gamma.margin <- function(param, y, x) {
  l <- length(param)
  eta <- x %*% param[-l]
  sum(stats::dgamma(y, shape = exp(param[l]), rate = exp(param[l])/exp(eta), log = TRUE))
}
loglik.cop <- function(param, u, copula) {
  if (class(copula)=="rotExplicitCopula"){
    copula@copula@parameters <- param
  } else {copula@parameters <- param}
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
         frankCopula  = { cop.param.condition <- ifelse((param.cop!=0), TRUE, FALSE) },
         rotExplicitCopula = {
           switch(class(copula@copula),
                  gumbelCopula = { cop.param.condition <- ifelse((param.cop>=1), TRUE, FALSE) },
                  joeCopula    = { cop.param.condition <- ifelse((param.cop>=1), TRUE, FALSE) },
                  claytonCopula= { cop.param.condition <- ifelse((param.cop>=-1 && param.cop!=0), TRUE, FALSE) }
           )
         }
  )
  if ( cop.param.condition ){
    u <- cbind(probtrans.gamma.margin(param1, y[,1], xmat[[1]]),
               probtrans.gamma.margin(param2, y[,2], xmat[[2]]))
    if (class(copula)=="rotExplicitCopula"){
      copula@copula@parameters <- param.cop
    } else {copula@parameters <- param.cop}
    loglik <- loglik.gamma.margin(param1, y[,1], xmat[[1]]) +
      loglik.gamma.margin(param2, y[,2], xmat[[2]]) +
      loglik.cop(param.cop, u, copula)
  } else { loglik <- NA }
  return(loglik)
}

thelogden.gamma <- function(theta, y, xmat, copula){
  l1 <- length(xmat[[1]]) + 1
  l2 <- length(xmat[[2]]) + 1
  param1 <- theta[1:l1]
  param2 <- theta[(l1 + 1):(l1 + l2)]
  param.cop <- theta[-(1:(l1 + l2))]
  switch(class(copula),
         tCopula      = { cop.param.condition <- ifelse((param.cop[1]<=1 && param.cop[1]>=-1), TRUE, FALSE) },
         normalCopula = { cop.param.condition <- ifelse((param.cop<=1 && param.cop>=-1), TRUE, FALSE) },
         gumbelCopula = { cop.param.condition <- ifelse((param.cop>=1), TRUE, FALSE) },
         joeCopula    = { cop.param.condition <- ifelse((param.cop>=1), TRUE, FALSE) },
         claytonCopula= { cop.param.condition <- ifelse((param.cop>=-1 && param.cop!=0), TRUE, FALSE) },
         frankCopula  = { cop.param.condition <- ifelse((param.cop!=0), TRUE, FALSE) },
         rotExplicitCopula = {
           switch(class(copula@copula),
                  gumbelCopula = { cop.param.condition <- ifelse((param.cop>=1), TRUE, FALSE) },
                  joeCopula    = { cop.param.condition <- ifelse((param.cop>=1), TRUE, FALSE) },
                  claytonCopula= { cop.param.condition <- ifelse((param.cop>=-1 && param.cop!=0), TRUE, FALSE) }
           )
         }
  )
  if ( cop.param.condition ){
    u <- cbind(probtrans.gamma.margin(param1, y[1], xmat[[1]]),
               probtrans.gamma.margin(param2, y[2], xmat[[2]]))
    if (class(copula)=="rotExplicitCopula"){
      copula@copula@parameters <- param.cop
    } else {copula@parameters <- param.cop}
    logden <- loglik.gamma.margin(param1, y[1], xmat[[1]]) +
              loglik.gamma.margin(param2, y[2], xmat[[2]]) +
              loglik.cop(param.cop, u, copula)
  } else { logden <- NA }
  return(logden)
}

Q.function <- function(theta, y, xmat, copula, z.mat){
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
         frankCopula  = { cop.param.condition <- ifelse((param.cop!=0), TRUE, FALSE) },
         rotExplicitCopula = {
           switch(class(copula@copula),
                  gumbelCopula = { cop.param.condition <- ifelse((param.cop>=1), TRUE, FALSE) },
                  joeCopula    = { cop.param.condition <- ifelse((param.cop>=1), TRUE, FALSE) },
                  claytonCopula= { cop.param.condition <- ifelse((param.cop>=-1 && param.cop!=0), TRUE, FALSE) }
                  )
         }
  )
  if ( cop.param.condition ){
    u <- cbind(probtrans.gamma.margin(param1, y[,1], xmat[[1]]),
               probtrans.gamma.margin(param2, y[,2], xmat[[2]]))
    if (class(copula)=="rotExplicitCopula"){
      copula@copula@parameters <- param.cop
    } else {copula@parameters <- param.cop}
    #first margin loglik
    eta1 <- xmat[[1]] %*% param1[-l1]
    loglik1 <- z.mat%*%stats::dgamma(y[,1], shape=exp(param1[l1]), rate=exp(param1[l1])/exp(eta1), log = TRUE)
    #second margin loglik
    eta2 <- xmat[[2]] %*% param2[-l2]
    loglik2 <- z.mat%*%stats::dgamma(y[,2], shape=exp(param2[l2]), rate=exp(param2[l2])/exp(eta2), log = TRUE)
    #copula loglik
    loglik3 <- z.mat%*%copula::dCopula(u, copula, log=TRUE )
    loglik <- loglik1+loglik2+loglik3
  } else { loglik <- NA }
  return(loglik)
}

