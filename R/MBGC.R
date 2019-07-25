#' Mixture of bivariate gamma distributions clustering
#'
#' Estimation using EM algorithm for mixture of bivariate gamma distributions
#'
#' @param modelName A character string indicating which model to be fitted. Need to be one of "CC", "CI", "IC".
#' @param y A vector of character strings indicating which variables in the data are treated as response or dependent variables.
#' @param G An integer specifying the numbers of mixture components.
#' @param gating Specifies the gating network in the MoE model, can be "C", "E" or a regression formula.
#' @param data A matrix or data frame of observations. Categorical variables are allowed as covariates.
#' @param maxit A parameter that controls the number of maximum iteration in the EM algorithm. The default is 100.
#' @param tol A parameter that controls the convergence tolerance in the EM algorithm. The default is 1e-5.
#' @param initialization Specifies initialization method for EM algorithm. The default is "\code{mclust}".
#' @param verbose A logical controlling whether estimations in each EM iteration are shown in the fitting process. The default is TRUE.
#'
#' @return An object of class \code{MBGC} providing the estimation results.
#'   The details of the output components are:
#'   \item{modelName}{A character string denoting the fitted expert model type.}
#'   \item{gating}{A character string denoting the fitted gating model type.}
#'   \item{alpha1}{The estimated alpha1 values.}
#'   \item{alpha2}{The estimated alpha2 values.}
#'   \item{alpha3}{The estimated alpha3 values.}
#'   \item{beta}{The estimated beta values.}
#'   \item{gating.coef}{The regression coefficients in the gating network, if a regression formula is provided in \code{gating}.}
#'   \item{pro}{A vector whose g-th component is the mixing proportion for the g-th component of the mixture model.}
#'   \item{z}{A matrix whose [i,g]-th entry is the probability that observation i in the data belongs to the g-th group.}
#'   \item{class}{The classification corresponding to z.}
#'   \item{G}{The number of mixture components.}
#'   \item{loglike}{The final estimated maximum log-likelihood value.}
#'   \item{ll}{The sequence of log-likelihood values in the EM algorithm fitting process.}
#'   \item{df}{Number of estimated parameters.}
#'   \item{AIC}{AIC values.}
#'   \item{BIC}{BIC values.}
#'   \item{iter}{Total iteration numbers.}
#'   \item{formula}{The formulas used in the regression.}
#'   \item{gating.model}{The final fitted regression model in gating network.}
#'   \item{y}{The input response data.}
#'   \item{n}{The number of observations in the data.}
#'   \item{Hessian}{The Hessian matrix at the estimated values}
#'   \item{call}{The matched call.}
#'   \item{trace}{All estimated coefficients and alpha, beta values in the EM algorithm.}
#'
#' @importFrom mclust Mclust mclustBIC
#'
#' @examples
#' \dontrun{
#' clust1 <- MBGC(modelName = "CC", y=c("y1","y2"),
#'                G=2, gating = "C", data=gatingsim, verbose=FALSE)
#' clust1
#' clust2 <- MBGC(modelName = "CI", y=c("y1","y2"),
#'                G=2, gating = "C", data=gatingsim)
#' clust2
#' clust3 <- MBGC(modelName = "IC", y=c("y1","y2"),
#'                G=2, gating = "C", data=gatingsim)
#' clust3
#' clust4 <- MBGC(modelName = "CC", y=c("y1","y2"),
#'                G=2, gating = ~w1+w2+w3, data=gatingsim)
#' clust4
#' clust5 <- MBGC(modelName = "CI", y=c("y1","y2"),
#'                G=2, gating = ~w1+w2+w3, data=gatingsim)
#' clust5
#' clust6 <- MBGC(modelName = "IC", y=c("y1","y2"),
#'                G=2, gating = ~w1+w2+w3, data=gatingsim)
#' clust6
#' }
#'
#' @importFrom stats uniroot Gamma formula coef optimHess
#' @importFrom mclust Mclust mclustBIC
#' @export

MBGC <- function(modelName = c("CC","CI","IC"),
                 y,
                 G,
                 gating,
                 data,
                 maxit = 100,
                 tol   = 1e-4,
                 initialization = "mclust",
                 verbose = TRUE){
  switch(modelName,
         CC = {
           res <- MBGC_CC(y=y, G=G, gating=gating, data=data,
                          maxit=maxit, tol=tol, initialization=initialization,
                          verbose=verbose)
         },
         CI = {
           res <- MBGC_CI(y=y, G=G, gating=gating, data= data,
                          maxit=maxit, tol=tol, initialization=initialization,
                          verbose=verbose)
         },
         IC = {
           res <- MBGC_IC(y=y, G=G, gating=gating, data=data,
                          maxit=maxit, tol=tol, initialization=initialization,
                          verbose=verbose)
         },
         stop("invalid model type name") )
  return(res)
}

#' @export

print.MBGC <- function (x, ...){
  modelfullname <- paste0(x$gating, x$modelName, collapse="")
  txt <- paste0("'", class(x)[1], "' model of type '", modelfullname, "' with G = ", x$G)
  cat(txt, "\n")
  cat("\n")
  cat("Available components:\n")
  print(names(x))
  invisible(x)
}

#' @rdname MBGC
#' @export
MBGC_CC <- function(y,
                    G,
                    gating, # "C", "E", formula
                    data,
                    maxit   = 200,
                    tol     = 1e-4,
                    initialization = "mclust",
                    verbose = TRUE){
  options(warn=-1)
  tempcall<-as.call(c(expression(MBGC.CC), list(y      = y,
                                                G      = G,
                                                gating = gating,
                                                data   = substitute(data),
                                                maxit  = maxit,
                                                tol    = tol,
                                                initialization = initialization,
                                                verbose= verbose)))

  namey1 <- y[1]
  namey2 <- y[2]
  y1     <- data[,names(data)==namey1]
  y2     <- data[,names(data)==namey2]
  n      <- length(y1)

  newdata <- data
  newdata <- newdata[ , names(newdata)!=namey1]
  newdata <- newdata[ , names(newdata)!=namey2]

  if (gating!="C" && gating!="E"){
    if (as.character(gating[2])=="."){gating.new <- formula( paste( "pzy", paste( names(newdata), "", collapse = "+", sep = ""), sep="~"))}
    else {gating.new <- formula(paste("pzy", as.character(gating)[2], sep="~"))}
  }

  #-----------------------------
  # starting values
  if (initialization=="mclust"){
    initial.mclust <- mclust::Mclust(G=G, data = cbind(y1,y2), verbose = FALSE)
    z.init      <- initial.mclust$z

    if (gating=="C"){
      p.z.current        <- initial.mclust$parameters$pro
      coef.p.current     <- NULL } else if (gating=="E"){
        p.z.current      <- rep(1/G,G)
        coef.p.current   <- NULL } else {
          newdata$pzy    <- z.init
          mp             <- nnet::multinom(formula = gating.new, data = newdata, trace=FALSE)  # reltol
          coef.p.current <- coef(mp)
          p.z.current    <- stats::predict(mp, type="probs") }

    index          <- initial.mclust$classification
    alpha1.current <- rep(0,G)
    alpha2.current <- rep(0,G)
    alpha3.current <- rep(0,G)
    beta.current   <- rep(0,G)

    for (gg in seq_len(G)){
      start.data <- cbind(y1,y2)[index==gg,]
      start.v1 <- MASS::fitdistr(start.data[,1]/mean(start.data[,1]), "gamma")
      start.v2 <- MASS::fitdistr(start.data[,2]/mean(start.data[,2]), "gamma")
      beta.current[gg] <- mean(c(start.v1$estimate[2]/mean(start.data[,1]),
                                 start.v2$estimate[2]/mean(start.data[,2])))
      alpha1.start <- start.v1$estimate[1]
      alpha2.start <- start.v2$estimate[1]
      alpha3.current[gg] <- min(alpha1.start, alpha2.start)/3
      alpha1.current[gg] <- alpha1.start - alpha3.current[gg]
      alpha2.current[gg] <- alpha2.start - alpha3.current[gg]
    }
  }

  j               <- 1
  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)

  while ( (loglike.diff > tol) && (j <= maxit) ) {
    #-------
    #E step
    #-------
    Expected.x3    <- matrix(0,nrow=n, ncol=G)
    Expected.x1    <- matrix(0,nrow=n, ncol=G)
    Expected.x2    <- matrix(0,nrow=n, ncol=G)
    Expected.logx1 <- matrix(0,nrow=n, ncol=G)
    Expected.logx2 <- matrix(0,nrow=n, ncol=G)
    Expected.logx3 <- matrix(0,nrow=n, ncol=G)
    den            <- matrix(0,nrow=n, ncol=G)

    for (i in seq_len(n)){
      for (g in seq_len(G)){
        expected.latent.res  <- expected_latent(c(y1[i],y2[i]),
                                                alpha=c(alpha1.current[g],
                                                        alpha2.current[g],
                                                        alpha3.current[g]),
                                                beta=beta.current[g])
        Expected.x3[i,g]     <- expected.latent.res$Expected.s
        Expected.x1[i,g]     <- (y1[i]-expected.latent.res$Expected.s)
        Expected.x2[i,g]     <- (y2[i]-expected.latent.res$Expected.s)
        Expected.logx3[i,g]  <- expected.latent.res$Expected.logs
        Expected.logx1[i,g]  <- expected.latent.res$Expected.logy1s
        Expected.logx2[i,g]  <- expected.latent.res$Expected.logy2s
        if (Expected.x1[i,g] <=0 || Expected.x2[i,g] <=0){
          Expected.x3[i,g]  <- min(y1[i], y2[i])/2
          Expected.x1[i,g]  <- y1[i]-Expected.x3[i,g]
          Expected.x2[i,g]  <- y2[i]-Expected.x3[i,g]
          warning("negative expected x1 or x2 at ", i, " obs at ", j, " iteration", "\n")
        }
        den[i,g]             <- expected.latent.res$denominator
      }
    }
    if (gating=="C" || gating=="E"){
      newden   <- sweep(den, 2, log(p.z.current), FUN="+") } else{
        newden <- den + log(p.z.current)  }
    denom   <- matrixStats::rowLogSumExps(newden)
    z.new   <- exp(sweep(newden, 1, denom, FUN="-"))

    loglike.new  <- sum(denom)
    loglike[j]   <- loglike.new
    loglike.diff <- abs(loglike.new-loglike.current)/(1+abs(loglike.new))
    if (verbose){
      cat(c("iter:", j, "; current loglike:", round(loglike.new,4), "; loglike change:", round(loglike.diff,5), "\n"))
    }

    #--------
    # M step
    #--------
    if (gating == "C"){
      p.z.new        <- colSums(z.new)/n
      coef.p.new     <- NULL
      mp             <- NULL } else if (gating == "E"){
        p.z.new      <- rep(1/G,G)
        coef.p.new   <- NULL
        mp           <- NULL } else {
          newdata$pzy  <- z.new
          mp         <- nnet::multinom(formula = gating.new, data = newdata, trace=FALSE)  # reltol
          coef.p.new <- coef(mp)
          p.z.new    <- stats::predict(mp, type="probs") }

    alpha1.new <- rep(0, G)
    alpha2.new <- rep(0, G)
    alpha3.new <- rep(0, G)
    beta.new   <- rep(0, G)
    for (g in seq_len(G)){
      beta.new[g] <- ((alpha1.current[g]+alpha2.current[g]+alpha3.current[g])*sum(z.new[,g])) / (z.new[,g]%*%(Expected.x1[,g]+Expected.x2[,g]+Expected.x3[,g]))
      alpha1.rootfun <- function(alpha1.var,wh){
        log(beta.new[wh])*sum(z.new[,wh]) - digamma(alpha1.var)*sum(z.new[,wh]) + z.new[,wh]%*%Expected.logx1[,g]
      }
      alpha1.new[g] <- stats::uniroot(alpha1.rootfun,wh=g,
                               lower=.Machine$double.eps, #sqrt(.Machine$double.eps),
                               upper=100000,
                               tol = sqrt(.Machine$double.xmin))$root
      alpha2.rootfun <- function(alpha2.var,wh){
        log(beta.new[wh])*sum(z.new[,wh]) - digamma(alpha2.var)*sum(z.new[,wh]) + z.new[,wh]%*%Expected.logx2[,g]
      }
      alpha2.new[g] <- stats::uniroot(alpha2.rootfun,wh=g,
                               lower=.Machine$double.eps,
                               upper=100000,
                               tol = sqrt(.Machine$double.xmin))$root
      alpha3.rootfun <- function(alpha3.var,wh){
        log(beta.new[wh])*sum(z.new[,wh]) - digamma(alpha3.var)*sum(z.new[,wh]) + z.new[,wh]%*%Expected.logx3[,g]
      }
      alpha3.new[g] <- stats::uniroot(alpha3.rootfun,wh=g,
                               lower=.Machine$double.eps,
                               upper=100000,
                               tol=sqrt(.Machine$double.xmin))$root
    }

    if (verbose){
      cat(c("alpha1:", round(alpha1.new,4), "\n"))
      cat(c("alpha2:", round(alpha2.new,4), "\n"))
      cat(c("alpha3:", round(alpha3.new,4), "\n"))
      cat(c("beta:", round(beta.new,4), "\n"))
      if (is.matrix(p.z.new)) { cat(c("p.z:", round(colMeans(p.z.new),4), "\n")) } else cat(c("p.z:", round(p.z.new,4), "\n"))
    }

    loglike.current<- loglike.new
    alpha1.current <- alpha1.new
    alpha2.current <- alpha2.new
    alpha3.current <- alpha3.new
    beta.current   <- beta.new
    p.z.current    <- p.z.new
    coef.p.current <- coef.p.new

    j              <- j+1
  }

  if (!all(loglike[1:(j-1)] == cummax(loglike[1:(j-1)]))){
    cat("Loglikelihood is not strictly monotonically increasing!", "\n")
  }

  Q.a.function <- function(a.var, wh, Expected.log.value, z.mat){
    q.a.res <- sum(z.mat[,wh]) * log(beta.current[wh]) * a.var -
      sum(z.mat[,wh]) * lgamma(a.var) +
      sum(z.mat[,wh] * Expected.log.value[,wh]) * a.var
    return(q.a.res)
  }
  Q.b.function <- function(b.var, wh, z.mat){
    q.b.res <- log(b.var) * sum(z.new[,wh]) * (alpha1.current[wh]+alpha2.current[wh]+alpha3.current[wh]) -
      b.var * (z.new[,wh]%*%(Expected.x1[,wh]+Expected.x2[,wh]+Expected.x3[,wh]))
    return(q.b.res)
  }
  hessian1 <- NULL
  hessian2 <- NULL
  hessian3 <- NULL
  hessian4 <- NULL
  for (g in c(1:G)){
    hessian1.temp  <- stats::optimHess(par   = alpha1.current[g],
                                       fn    = Q.a.function,
                                       wh    = g,
                                       Expected.log.value = Expected.logx1,
                                       z.mat = z.new)
    hessian1       <- append(hessian1, list(hessian1.temp))
    hessian2.temp  <- stats::optimHess(par   = alpha2.current[g],
                                       fn    = Q.a.function,
                                       wh    = g,
                                       Expected.log.value = Expected.logx2,
                                       z.mat = z.new)
    hessian2       <- append(hessian2, list(hessian2.temp))
    hessian3.temp  <- stats::optimHess(par   = alpha3.current[g],
                                       fn    = Q.a.function,
                                       wh    = g,
                                       Expected.log.value = Expected.logx3,
                                       z.mat = z.new)
    hessian3       <- append(hessian3, list(hessian3.temp))
    hessian4.temp  <- stats::optimHess(par   = beta.current[g],
                                       fn    = Q.b.function,
                                       wh    = g,
                                       z.mat = z.new)
    hessian4       <- append(hessian4, list(hessian4.temp))
  }
  all.hessian <- c("alpha1,g="=hessian1, "alpha2,g="=hessian2, "alpha3,g="=hessian3, "beta,g="=hessian4)
  if (!all(sapply(all.hessian, matrixcalc::is.negative.definite))){
    cat("Hessian matrix is not negative-definite.", "\n")
  }

  #	calculation of BIC and AIC for bivpoisson model
  if (gating=="C"){
    noparams <- 5*G-1
  } else if (gating=="E"){
    noparams <- 4*G
  } else {
    noparams <- 4*G+mp$rank
  }
  AIC<- -2*loglike[j-1] + noparams * 2
  BIC<- -2*loglike[j-1] + noparams * log(n)
  classification <- apply(z.new, 1, which.max)

  if (gating == "C") {
    newgating <- gating.formula <- "C" } else if (gating == "E") {
      newgating <- gating.formula <- "E"} else {
        gating.formula <- formula(paste("", as.character(gating.new)[3], sep="~"))
        newgating <- "V"}

  result<-list(modelName    = "CC",
               gating       = newgating,
               alpha1       = alpha1.current,
               alpha2       = alpha2.current,
               alpha3       = alpha3.current,
               beta         = beta.current,
               gating.coef  = list(coef.p.current),
               pro          = p.z.current,
               z            = z.new,
               class        = classification,
               G            = G,
               loglike      = loglike.current,
               ll           = loglike[1:(j-1)],
               df           = noparams,
               AIC          = AIC,
               BIC          = BIC,
               n            = n,
               y            = cbind(y1, y2),
               formula      = gating.formula,
               gating.model = mp,
               Hessian      = all.hessian,
               call         = tempcall,
               iter         = (j-1))
  options(warn=0)
  structure(result, class = 'MBGC')
}

#' @rdname MBGC
#' @export
MBGC_CI <- function(y,
                    G,
                    gating, # "C", "E", formula
                    data,
                    maxit   = 200,
                    tol     = 1e-4,
                    initialization = "mclust",
                    verbose = TRUE){
  options(warn=-1)
  tempcall<-as.call(c(expression(MBGC.CI),list(y      = y,
                                               G      = G,
                                               gating = gating,
                                               data   = substitute(data),
                                               maxit  = maxit,
                                               tol    = tol,
                                               initialization = initialization,
                                               verbose= verbose)))
  namey1 <- y[1]
  namey2 <- y[2]
  y1     <- data[,names(data)==namey1]
  y2     <- data[,names(data)==namey2]
  n      <- length(y1)

  newdata <- data
  newdata <- newdata[ , names(newdata)!=namey1]
  newdata <- newdata[ , names(newdata)!=namey2]

  if (gating!="C" && gating!="E"){
    if (as.character(gating[2])=="."){gating.new <- formula( paste( "pzy", paste( names(newdata), "", collapse = "+", sep = ""), sep="~"))}
    else {gating.new <- formula(paste("pzy", as.character(gating)[2], sep="~"))}
  }

  #-----------------------------
  # starting values
  if (initialization=="mclust"){
    initial.mclust <- mclust::Mclust(G=G, data = cbind(y1,y2), verbose = FALSE)
    z.init      <- initial.mclust$z

    if (gating=="C"){
      p.z.current        <- initial.mclust$parameters$pro
      coef.p.current     <- NULL } else if (gating=="E"){
        p.z.current      <- rep(1/G,G)
        coef.p.current   <- NULL } else {
          newdata$pzy    <- z.init
          mp             <- nnet::multinom(formula = gating.new, data = newdata, trace=FALSE)  # reltol
          coef.p.current <- coef(mp)
          p.z.current    <- stats::predict(mp, type="probs") }

    index          <- initial.mclust$classification
    alpha1.current <- rep(0,G)
    alpha2.current <- rep(0,G)
    alpha3.current <- rep(0,G)
    beta.temp   <- rep(0,G)

    for (gg in seq_len(G)){
      start.data <- cbind(y1,y2)[index==gg,]
      start.v1 <- MASS::fitdistr(start.data[,1]/mean(start.data[,1]), "gamma")
      start.v2 <- MASS::fitdistr(start.data[,2]/mean(start.data[,2]), "gamma")
      beta.temp[gg] <- mean(c(start.v1$estimate[2]/mean(start.data[,1]),
                              start.v2$estimate[2]/mean(start.data[,2])))
      alpha1.start <- start.v1$estimate[1]
      alpha2.start <- start.v2$estimate[1]
      alpha3.current[gg] <- min(alpha1.start, alpha2.start)/3
      alpha1.current[gg] <- alpha1.start - alpha3.current[gg]
      alpha2.current[gg] <- alpha2.start - alpha3.current[gg]
    }
    beta.current <- mean(beta.temp)
  }

  j               <- 1
  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)

  while ( (loglike.diff > tol) && (j <= maxit) ) {
    #-------------------
    #E step
    #-------------------
    Expected.x3    <- matrix(0,nrow=n, ncol=G)
    Expected.x1    <- matrix(0,nrow=n, ncol=G)
    Expected.x2    <- matrix(0,nrow=n, ncol=G)
    Expected.logx1 <- matrix(0,nrow=n, ncol=G)
    Expected.logx2 <- matrix(0,nrow=n, ncol=G)
    Expected.logx3 <- matrix(0,nrow=n, ncol=G)
    den            <- matrix(0,nrow=n, ncol=G)
    for (i in seq_len(n)){
      for (g in seq_len(G)){
        expected.latent.res  <- expected_latent(c(y1[i],y2[i]),
                                                alpha=c(alpha1.current[g],
                                                        alpha2.current[g],
                                                        alpha3.current[g]),
                                                beta=beta.current)
        Expected.x3[i,g]     <- expected.latent.res$Expected.s
        Expected.x1[i,g]     <- (y1[i]-expected.latent.res$Expected.s)
        Expected.x2[i,g]     <- (y2[i]-expected.latent.res$Expected.s)
        Expected.logx3[i,g]  <- expected.latent.res$Expected.logs
        Expected.logx1[i,g]  <- expected.latent.res$Expected.logy1s
        Expected.logx2[i,g]  <- expected.latent.res$Expected.logy2s
        if (Expected.x1[i,g] <=0 || Expected.x2[i,g] <=0){
          Expected.x3[i,g]  <- min(y1[i], y2[i])/2
          Expected.x1[i,g]  <- y1[i]-Expected.x3[i,g]
          Expected.x2[i,g]  <- y2[i]-Expected.x3[i,g]
          warning("negative expected x1 or x2 at ", i, " obs at ", j, " iteration", "\n")
        }
        den[i,g]             <- expected.latent.res$denominator
      }
    }
    if (gating=="C" || gating=="E"){
      newden   <- sweep(den, 2, log(p.z.current), FUN="+") } else{
        newden <- den + log(p.z.current)  }
    denom   <- matrixStats::rowLogSumExps(newden)
    z.new   <- exp(sweep(newden, 1, denom, FUN="-"))

    loglike.new    <- sum(denom)
    loglike[j]     <- loglike.new
    loglike.diff <- abs(loglike.new-loglike.current)/(1+abs(loglike.new))
    if (verbose){
      cat(c("iter:", j, "; current loglike:", round(loglike.new,4), "; loglike change:", round(loglike.diff,5), "\n"))
    }

    #--------
    # M step
    #--------
    if (gating == "C"){
      p.z.new        <- colSums(z.new)/n
      coef.p.new     <- NULL
      mp             <- NULL } else if (gating == "E"){
        p.z.new      <- rep(1/G,G)
        coef.p.new   <- NULL
        mp           <- NULL } else {
          newdata$pzy  <- z.new
          mp         <- nnet::multinom(formula = gating.new, data = newdata, trace=FALSE)  # reltol
          coef.p.new <- coef(mp)
          p.z.new    <- stats::predict(mp, type="probs") }

    beta.new   <- (alpha1.current+alpha2.current+alpha3.current)%*%colSums(z.new) / sum(z.new*(Expected.x1+Expected.x2+Expected.x3))
    alpha1.new <- rep(0, G)
    alpha2.new <- rep(0, G)
    alpha3.new <- rep(0, G)
    for (g in seq_len(G)){
      alpha1.rootfun <- function(alpha1.var,wh){
        log(beta.new)*sum(z.new[,wh]) - digamma(alpha1.var)*sum(z.new[,wh]) + z.new[,wh]%*%Expected.logx1[,wh]
      }
      alpha1.new[g] <- stats::uniroot(alpha1.rootfun,wh=g,
                               lower=.Machine$double.eps,
                               upper=100000,
                               tol = sqrt(.Machine$double.xmin))$root
      alpha2.rootfun <- function(alpha2.var,wh){
        log(beta.new)*sum(z.new[,wh]) - digamma(alpha2.var)*sum(z.new[,wh]) + z.new[,wh]%*%Expected.logx2[,wh]
      }
      alpha2.new[g] <- stats::uniroot(alpha2.rootfun,wh=g,
                               lower=.Machine$double.eps,
                               upper=100000,
                               tol = sqrt(.Machine$double.xmin))$root
      alpha3.rootfun <- function(alpha3.var,wh){
        log(beta.new)*sum(z.new[,wh]) - digamma(alpha3.var)*sum(z.new[,wh]) + z.new[,wh]%*%Expected.logx3[,wh]
      }
      alpha3.new[g] <- stats::uniroot(alpha3.rootfun,wh=g,
                               lower=.Machine$double.eps,
                               upper=100000,
                               tol=sqrt(.Machine$double.xmin))$root
    }

    if (verbose){
      cat(c("alpha1:", round(alpha1.new,4), "\n"))
      cat(c("alpha2:", round(alpha2.new,4), "\n"))
      cat(c("alpha3:", round(alpha3.new,4), "\n"))
      cat(c("beta:", round(beta.new,4), "\n"))
      if (is.matrix(p.z.new)) { cat(c("p.z:", round(colMeans(p.z.new),4), "\n")) } else cat(c("p.z:", round(p.z.new,4), "\n"))
    }

    loglike.current<- loglike.new
    alpha1.current <- alpha1.new
    alpha2.current <- alpha2.new
    alpha3.current <- alpha3.new
    beta.current   <- beta.new
    p.z.current    <- p.z.new
    coef.p.current <- coef.p.new
    j              <- j+1
  }

  if (!all(loglike[1:(j-1)] == cummax(loglike[1:(j-1)]))){
    cat("Loglikelihood is not strictly monotonically increasing!", "\n")
  }

  Q.a.function <- function(a.var, wh, Expected.log.value, z.mat){
    q.a.res <- sum(z.mat[,wh]) * log(beta.current) * a.var -
      sum(z.mat[,wh]) * lgamma(a.var) +
      sum(z.mat[,wh] * Expected.log.value[,wh]) * a.var
    return(q.a.res)
  }
  Q.b.function <- function(b.var, z.mat){
    q.b.res <- log(b.var) * colSums(z.new)%*%(alpha1.current+alpha2.current+alpha3.current) -
      b.var * sum(z.new*(Expected.x1+Expected.x2+Expected.x3))
    return(q.b.res)
  }
  hessian1 <- NULL
  hessian2 <- NULL
  hessian3 <- NULL
  for (g in c(1:G)){
    hessian1.temp  <- stats::optimHess(par   = alpha1.current[g],
                                       fn    = Q.a.function,
                                       wh    = g,
                                       Expected.log.value = Expected.logx1,
                                       z.mat = z.new)
    hessian1       <- append(hessian1, list(hessian1.temp))
    hessian2.temp  <- stats::optimHess(par   = alpha2.current[g],
                                       fn    = Q.a.function,
                                       wh    = g,
                                       Expected.log.value = Expected.logx2,
                                       z.mat = z.new)
    hessian2       <- append(hessian2, list(hessian2.temp))
    hessian3.temp  <- stats::optimHess(par   = alpha3.current[g],
                                       fn    = Q.a.function,
                                       wh    = g,
                                       Expected.log.value = Expected.logx3,
                                       z.mat = z.new)
    hessian3       <- append(hessian3, list(hessian3.temp))
  }
  hessian4 <- stats::optimHess(par   = beta.current,
                               fn    = Q.b.function,
                               z.mat = z.new)
  all.hessian <- c("alpha1,g="=hessian1, "alpha2,g="=hessian2, "alpha3,g="=hessian3, "beta"=list(hessian4))
  if (!all(sapply(all.hessian, matrixcalc::is.negative.definite))){
    cat("Hessian matrix is not negative-definite.", "\n")
  }

  #	calculation of BIC and AIC for bivpoisson model
  if (gating=="C"){
    noparams <- 4*G
  } else if (gating=="E"){
    noparams <- 3*G+1
  } else {
    noparams <- 3*G+1+mp$rank
  }
  AIC<- -2*loglike[j-1] + noparams * 2
  BIC<- -2*loglike[j-1] + noparams * log(n)
  classification <- apply(z.new, 1, which.max)

  if (gating == "C") {
    newgating <- gating.formula <- "C" } else if (gating == "E") {
      newgating <- gating.formula <- "E"} else {
        gating.formula <- formula(paste("", as.character(gating.new)[3], sep="~"))
        newgating <- "V" }

  result<-list(modelName    = "CI",
               gating       = newgating,
               alpha1       = alpha1.current,
               alpha2       = alpha2.current,
               alpha3       = alpha3.current,
               beta         = beta.current,
               gating.coef  = list(coef.p.current),
               pro          = p.z.current,
               z            = z.new,
               G            = G,
               class        = classification,
               loglike      = loglike.current,
               ll           = loglike[1:(j-1)],
               df           = noparams,
               AIC          = AIC,
               BIC          = BIC,
               n            = n,
               y            = cbind(y1, y2),
               formula      = gating.formula,
               gating.model = mp,
               Hessian      = all.hessian,
               call         = tempcall,
               iter         = (j-1))
  options(warn=0)
  structure(result, class = 'MBGC')
}

#' @rdname MBGC
#' @export
MBGC_IC <- function(y,
                    G,
                    gating, # "C", "E", formula
                    data,
                    maxit   = 200,
                    tol     = 1e-4,
                    initialization = "mclust",
                    verbose = TRUE){
  options(warn=-1)
  tempcall<-as.call(c(expression(MBGC), list(y      = y,
                                             G      = G,
                                             gating = gating,
                                             data   = substitute(data),
                                             maxit  = maxit,
                                             tol    = tol,
                                             initialization = initialization,
                                             verbose= verbose)))
  namey1 <- y[1]
  namey2 <- y[2]
  y1     <- data[,names(data)==namey1]
  y2     <- data[,names(data)==namey2]
  n      <- length(y1)

  newdata <- data
  newdata <- newdata[ , names(newdata)!=namey1]
  newdata <- newdata[ , names(newdata)!=namey2]

  if (gating!="C" && gating!="E"){
    if (as.character(gating[2])=="."){gating.new <- formula( paste( "pzy", paste( names(newdata), "", collapse = "+", sep = ""), sep="~"))}
    else {gating.new <- formula(paste("pzy", as.character(gating)[2], sep="~"))}
  }

  #-----------------------------
  # starting values
  if (initialization=="mclust"){
    initial.mclust <- Mclust(G=G, data = cbind(y1,y2), verbose = FALSE)
    z.init      <- initial.mclust$z

    if (gating=="C"){
      p.z.current        <- initial.mclust$parameters$pro
      coef.p.current     <- NULL } else if (gating=="E"){
        p.z.current      <- rep(1/G,G)
        coef.p.current   <- NULL } else {
          newdata$pzy    <- z.init
          mp             <- nnet::multinom(formula = gating.new, data = newdata, trace=FALSE)  # reltol
          coef.p.current <- coef(mp)
          p.z.current    <- stats::predict(mp, type="probs") }

    index          <- initial.mclust$classification
    alpha1.temp <- rep(0,G)
    alpha2.temp <- rep(0,G)
    alpha3.temp <- rep(0,G)
    beta.current   <- rep(0,G)

    for (gg in seq_len(G)){
      start.data <- cbind(y1,y2)[index==gg,]
      start.v1 <- MASS::fitdistr(start.data[,1]/mean(start.data[,1]), "gamma")
      start.v2 <- MASS::fitdistr(start.data[,2]/mean(start.data[,2]), "gamma")
      beta.current[gg] <- mean(c(start.v1$estimate[2]/mean(start.data[,1]),
                                 start.v2$estimate[2]/mean(start.data[,2])))
      alpha1.start <- start.v1$estimate[1]
      alpha2.start <- start.v2$estimate[1]
      alpha3.temp[gg] <- min(alpha1.start, alpha2.start)/3
      alpha1.temp[gg] <- alpha1.start - alpha3.temp[gg]
      alpha2.temp[gg] <- alpha2.start - alpha3.temp[gg]
    }
    alpha1.current <- mean(alpha1.temp)
    alpha2.current <- mean(alpha2.temp)
    alpha3.current <- mean(alpha3.temp)
  }

  j               <- 1
  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)

  while ( (loglike.diff > tol) && (j <= maxit) ) {
    #-------
    #E step
    #-------
    Expected.x3    <- matrix(0,nrow=n, ncol=G)
    Expected.x1    <- matrix(0,nrow=n, ncol=G)
    Expected.x2    <- matrix(0,nrow=n, ncol=G)
    Expected.logx1 <- matrix(0,nrow=n, ncol=G)
    Expected.logx2 <- matrix(0,nrow=n, ncol=G)
    Expected.logx3 <- matrix(0,nrow=n, ncol=G)
    den            <- matrix(0,nrow=n, ncol=G)

    for (i in seq_len(n)){
      for (g in seq_len(G)){
        expected.latent.res  <- expected_latent(c(y1[i],y2[i]),
                                                alpha=c(alpha1.current,
                                                        alpha2.current,
                                                        alpha3.current),
                                                beta=beta.current[g])
        Expected.x3[i,g]     <- expected.latent.res$Expected.s
        Expected.x1[i,g]     <- (y1[i]-expected.latent.res$Expected.s)
        Expected.x2[i,g]     <- (y2[i]-expected.latent.res$Expected.s)
        Expected.logx3[i,g]  <- expected.latent.res$Expected.logs
        Expected.logx1[i,g]  <- expected.latent.res$Expected.logy1s
        Expected.logx2[i,g]  <- expected.latent.res$Expected.logy2s
        if (Expected.x1[i,g] <=0 || Expected.x2[i,g] <=0){
          Expected.x3[i,g]  <- min(y1[i], y2[i])/2
          Expected.x1[i,g]  <- y1[i]-Expected.x3[i,g]
          Expected.x2[i,g]  <- y2[i]-Expected.x3[i,g]
          warning("negative expected x1 or x2 at ", i, " obs at ", j, " iteration", "\n")
        }
        den[i,g]             <- expected.latent.res$denominator
      }
    }
    if (gating=="C" || gating=="E"){
      newden   <- sweep(den, 2, log(p.z.current), FUN="+") } else{
        newden <- den + log(p.z.current)  }
    denom   <- matrixStats::rowLogSumExps(newden)
    z.new   <- exp(sweep(newden, 1, denom, FUN="-"))

    loglike.new    <- sum(denom)
    loglike[j]     <- loglike.new
    loglike.diff <- abs(loglike.new-loglike.current)/(1+abs(loglike.new))
    if (verbose){
      cat(c("iter:", j, "; current loglike:", round(loglike.new,4), "; loglike change:", round(loglike.diff,5), "\n"))
    }

    #--------
    # M step
    #--------
    if (gating == "C"){
      p.z.new        <- colSums(z.new)/n
      coef.p.new     <- NULL
      mp             <- NULL } else if (gating == "E"){
        p.z.new      <- rep(1/G,G)
        coef.p.new   <- NULL
        mp           <- NULL } else {
          newdata$pzy  <- z.new
          mp         <- nnet::multinom(formula = gating.new, data = newdata, trace=FALSE)  # reltol
          coef.p.new <- coef(mp)
          p.z.new    <- stats::predict(mp, type="probs") }

    beta.new <- ((alpha1.current+alpha2.current+alpha3.current)*colSums(z.new)) / colSums(z.new*(Expected.x1+Expected.x2+Expected.x3))
    alpha1.rootfun <- function(alpha1.var){
      sum(log(beta.new)*colSums(z.new)) - digamma(alpha1.var)*sum(z.new) + sum(z.new*Expected.logx1)
    }
    alpha1.new <- stats::uniroot(alpha1.rootfun,
                          lower=.Machine$double.eps,
                          upper=100000,
                          tol = sqrt(.Machine$double.xmin))$root
    alpha2.rootfun <- function(alpha2.var){
      sum(log(beta.new)*colSums(z.new)) - digamma(alpha2.var)*sum(z.new) + sum(z.new*Expected.logx2)
    }
    alpha2.new <- stats::uniroot(alpha2.rootfun,
                          lower=.Machine$double.eps,
                          upper=100000,
                          tol = sqrt(.Machine$double.xmin))$root
    alpha3.rootfun <- function(alpha3.var){
      sum(log(beta.new)*colSums(z.new)) - digamma(alpha3.var)*sum(z.new) + sum(z.new*Expected.logx3)
    }
    alpha3.new <- stats::uniroot(alpha3.rootfun,
                          lower=.Machine$double.eps,
                          upper=100000,
                          tol=sqrt(.Machine$double.xmin))$root

    if (verbose){
      cat(c("alpha1:", round(alpha1.new,4), "\n"))
      cat(c("alpha2:", round(alpha2.new,4), "\n"))
      cat(c("alpha3:", round(alpha3.new,4), "\n"))
      cat(c("beta:", round(beta.new,4), "\n"))
      if (is.matrix(p.z.new)) { cat(c("p.z:", round(colMeans(p.z.new),4), "\n")) } else cat(c("p.z:", round(p.z.new,4), "\n"))
    }

    loglike.current<- loglike.new
    alpha1.current <- alpha1.new
    alpha2.current <- alpha2.new
    alpha3.current <- alpha3.new
    beta.current   <- beta.new
    p.z.current    <- p.z.new
    coef.p.current <- coef.p.new
    j              <- j+1
  }

  if (!all(loglike[1:(j-1)] == cummax(loglike[1:(j-1)]))){
    cat("Loglikelihood is not strictly monotonically increasing!", "\n")
  }

  Q.a.function <- function(a.var, Expected.log.value){
    q.a.res <- sum(log(beta.new)*colSums(z.new)) * a.var -
      sum(z.new) * lgamma(a.var) +
      sum(z.new * Expected.log.value) * a.var
    return(q.a.res)
  }
  Q.b.function <- function(b.var, wh){
    q.b.res <- log(b.var) * sum(z.new[,wh]) * (alpha1.current+alpha2.current+alpha3.current) -
      b.var * (z.new[,wh]%*%(Expected.x1[,wh]+Expected.x2[,wh]+Expected.x3[,wh]))
    return(q.b.res)
  }
  hessian1 <- stats::optimHess(par   = alpha1.current,
                               fn    = Q.a.function,
                               Expected.log.value = Expected.logx1)
  hessian2 <- stats::optimHess(par   = alpha2.current,
                               fn    = Q.a.function,
                               Expected.log.value = Expected.logx2)
  hessian3 <- stats::optimHess(par   = alpha3.current,
                               fn    = Q.a.function,
                               Expected.log.value = Expected.logx3)
  hessian4 <- NULL
  for (g in c(1:G)){
    hessian4.temp  <- stats::optimHess(par   = beta.current[g],
                                       fn    = Q.b.function,
                                       wh    = g)
    hessian4       <- append(hessian4, list(hessian4.temp))
  }
  all.hessian <- c("alpha1"=list(hessian1), "alpha2"=list(hessian2), "alpha3"=list(hessian3), "beta,g="=hessian4)
  if (!all(sapply(all.hessian, matrixcalc::is.negative.definite))){
    cat("Hessian matrix is not negative-definite.", "\n")
  }

  #	calculation of BIC and AIC for bivpoisson model
  if (gating=="C"){
    noparams <- 2*G+2
  } else if (gating=="E"){
    noparams <- G+3
  } else {
    noparams <- G+3+mp$rank
  }
  AIC<- -2*loglike[j-1] + noparams * 2
  BIC<- -2*loglike[j-1] + noparams * log(n)
  classification <- apply(z.new, 1, which.max)

  if (gating == "C") {
    newgating <- gating.formula <- "C" } else if (gating == "E") {
      newgating <- gating.formula <- "E"} else {
        gating.formula <- formula(paste("", as.character(gating.new)[3], sep="~"))
        newgating <- "V"}

  result<-list(modelName    = "IC",
               gating       = newgating,
               alpha1       = alpha1.current,
               alpha2       = alpha2.current,
               alpha3       = alpha3.current,
               beta         = beta.current,
               gating.coef  = list(coef.p.current),
               pro          = p.z.current,
               z            = z.new,
               class        = classification,
               G            = G,
               loglike      = loglike.current,
               ll           = loglike[1:(j-1)],
               df           = noparams,
               AIC          = AIC,
               BIC          = BIC,
               n            = n,
               y            = cbind(y1, y2),
               formula      = gating.formula,
               gating.model = mp,
               Hessian      = all.hessian,
               call         = tempcall,
               iter         = (j-1))
  options(warn=0)
  structure(result, class = 'MBGC')
}



