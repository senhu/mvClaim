#' Tail Dependence Coefficient
#'
#' This function returns the tail dependence coefficients that estimated nonparametrically using empirical copulas.
#'
#' @param X A data matrix
#' @param lower Logical; if \code{TRUE} the lower TDC is computed, and upper TDC otherwise.
#' @param method Specifies which method to use, integer values from 1 to 4.
#' @param plot Logical; if \code{TRUE} the function also plots the corresponding plot during estimation. Only for method 3 and 4.
#' @param control A list of control parameters for method 4.
#'
#' @return The estimated tail dependence coefficient value
#'
#' @examples
#' data(simdat.mcgr)
#' TDC(simdat.mcgr[,1:2], method = 1, lower=FALSE)
#' TDC(simdat.mcgr[,1:2], method = 2, lower=FALSE)
#' TDC(simdat.mcgr[,1:2], method = 3, lower=FALSE)
#'
#' @importFrom graphics lines abline
#' @export TDC

TDC <- function(X, lower, method, plot, control){
  switch(method,
         "1" = {
           tdc.empcop.1(X, lower)
         },
         "2" = {
           tdc.empcop.2(X, lower)
         },
         "3" = {
           tdc.empcop.3(X, lower, plot = FALSE)
         },
         "4" = {
           tdc.empcop.4(X, lower, plot = TRUE,
                        control = list(bootstrap = FALSE,
                                       bootstrap.num = 100,
                                       smooth = TRUE,
                                       cutoff = 0.9,
                                       smoothing.span = 0.1))
         }, stop("invalid input method of TDC estimation"))
}


tdc.empcop.1 <- function(X, lower = FALSE){
  n <- nrow(X)
  t <- sqrt(n)
  if (!lower) {
    c.approx <- copula::C.n(matrix(c(1-t/n, 1-t/n),nrow=1), X)
    utdc <- (1 - 2*(1-t/n) + c.approx)/(t/n)
    return(utdc)
  } else {
    c.approx <- copula::C.n(matrix(c(t/n, t/n),nrow=1), X, smoothing = "none")
    ltdc <- c.approx/(t/n)
    return(ltdc)
  }
}

tdc.empcop.2 <- function(X, lower){
  n <- nrow(X)
  t <- sqrt(n)
  if (!lower) {
    c.approx <- copula::C.n(matrix(c(1-t/n, 1-t/n),nrow=1), X, smoothing = "none")
    utdc <- (2 - log(c.approx)/log(1-t/n))
    return(utdc)
  } else {
    c.approx <- copula::C.n(matrix(c(t/n, t/n),nrow=1), X, smoothing = "none")
    ltdc <- 2 - log(1 - 2*(t/n) + c.approx)/log(1-t/n)
    return(ltdc)
  }
}

tdc.empcop.3 <- function(X, lower, plot=FALSE){
  n <- nrow(X)
  k <- round(sqrt(n))
  if (!lower){
    Uemp <- function(i) {1 - 2*(1-i/n) + copula::C.n(matrix(c(1-i/n, 1-i/n),nrow=1), X, smoothing = "none")}
    emp.lambda <- Vectorize(Uemp)(k:1)
    xx <- (k:1)/n
    lmod <- stats::lm(emp.lambda~(xx-1))
    if (plot) {
      plot(xx, emp.lambda, type="l", ylab="upper tail dependence")
      lines(xx, lmod$fitted.values, col=2)
    }
    return(sum(coef(lmod)))
  } else {
    Lemp <- function(i) {copula::C.n(matrix(c(i/n, i/n),nrow=1), X, smoothing = "none")}
    emp.lambda <- Vectorize(Lemp)(1:k)
    xx <- (1:k)/n
    lmod <- stats::lm(emp.lambda~(xx-1))
    if (plot) {
      plot(xx, emp.lambda, type="l", ylab="lower tail dependence")
      lines(xx, lmod$fitted.values, col=2)
    }
    return(sum(coef(lmod)))
  }
}

tdc.empcop.4 <- function(X, lower, plot = TRUE,
                         control = list(bootstrap = FALSE,
                                        bootstrap.num = 100,
                                        smooth = TRUE,
                                        cutoff = 0.9,
                                        smoothing.span = 0.1)){
  if (!lower){
    n <- nrow(X)
    if (!control$bootstrap){
      U=X[order(X[,1]),1]
      V=X[order(X[,2]),2]
      Emp <- function(i, j) {sum(X[,1]<=U[i] & X[,2]<=V[j])/n}
      Uemp <- function(i) (1-(2*i/n)+Emp(i, i))/(1-(i/n))
      emp.lamb <- Vectorize(Uemp)(1:(n-1))
      if (!control$smooth){
        wm <- which.min(emp.lamb[1:((n-1)*control$cutoff)])
        m  <- min(emp.lamb[1:((n-1)*control$cutoff)])
      }
      if (control$smooth) {
        con <- 1:(n-1)
        emp.lamb.ss <- stats::predict(stats::loess(emp.lamb ~ con, span = control$smoothing.span))
        wm <- which.min(emp.lamb.ss[1:((n-1)*control$cutoff)])
        m <- min(emp.lamb.ss[1:((n-1)*control$cutoff)])
      }
      if (plot){
        plot(emp.lamb, ylim=c(0,1), type="l", ylab="upper tail dependence")
        if (control$smooth) lines(emp.lamb.ss, col=2)
        abline(h=m, v=wm)
      }
    }
    if (control$bootstrap){
      boot.res <- NULL
      for (k in 1:control$bootstrap.num){
        Z <- X[sample(1:n, n, replace=TRUE),]
        U.z=Z[order(Z[,1]),1]
        V.z=Z[order(Z[,2]),2]
        Emp.boot <- function(i, j) {sum(Z[,1]<=U.z[i] & Z[,2]<=V.z[j])/n}
        Uemp.boot <- function(i) (1-(2*i/n)+Emp.boot(i, i))/(1-(i/n))
        boot.res <- rbind(boot.res, Vectorize(Uemp.boot)(1:(n-1)))
      }
      lm <- apply(boot.res, 2, mean)
      lq <- apply(boot.res[,1:(n-1)], 2, stats::quantile, .025)
      uq <- apply(boot.res[,1:(n-1)], 2, stats::quantile, .975)
      if (!control$smooth){
        wm <- which.min(lm[1:((n-1)*control$cutoff)])
        m  <- min(lm[1:((n-1)*control$cutoff)])
      }
      if (control$smooth) {
        con <- 1:(n-1)
        emp.lamb.ss <- stats::predict(stats::loess(lm ~ con, span = control$smoothing.span))
        wm <- which.min(emp.lamb.ss[1:((n-1)*control$cutoff)])
        m <- min(emp.lamb.ss[1:((n-1)*control$cutoff)])
      }
      if (plot){
        plot(boot.res[1,], ylim=c(0,1), type="l", ylab="upper tail dependence")
        for (k in 2:100) lines(boot.res[k,])
        lines(lm, col=2, lwd=2)
        lines(lq, col=2)
        lines(uq, col=2)
        if (control$smooth) lines(emp.lamb.ss, col=4, lwd=2)
        abline(h=m, v=wm)
      }
    }
    return(c(m, wm))
  } else {
    n <- dim(X)[1]
    if (!control$bootstrap){
      U=X[order(X[,1]),1]
      V=X[order(X[,2]),2]
      Emp <- function(i, j) {sum(X[,1]<=U[i] & X[,2]<=V[j])/n}
      Lemp <- function(i) Emp(i, i)/(i/n)
      emp.lamb <- Vectorize(Lemp)(1:(n-1))
      if (!control$smooth){
        wm <- which.min(emp.lamb[((n-1)*(1-control$cutoff)):(n-1)]) + (n-1)*(1-control$cutoff)
        m  <- min(emp.lamb[((n-1)*(1-control$cutoff)):(n-1)])
      }
      if (control$smooth) {
        con <- 1:(n-1)
        emp.lamb.ss <- stats::predict(stats::loess(emp.lamb ~ con, span = control$smoothing.span))
        wm <- which.min(emp.lamb.ss[((n-1)*(1-control$cutoff)):(n-1)]) + (n-1)*(1-control$cutoff)
        m <- min(emp.lamb.ss[((n-1)*(1-control$cutoff)):(n-1)])
      }
      if (plot){
        plot(emp.lamb, ylim=c(0,1), type="l", ylab="upper tail dependence")
        if (control$smooth) lines(emp.lamb.ss, col=2)
        abline(h=m, v=wm)
      }
    }
    if (control$bootstrap){
      boot.res <- NULL
      for (k in 1:control$bootstrap.num){
        Z <- X[sample(1:n, n, replace=TRUE),]
        U.z=Z[order(Z[,1]),1]
        V.z=Z[order(Z[,2]),2]
        Emp.boot <- function(i, j) {sum(Z[,1]<=U.z[i] & Z[,2]<=V.z[j])/n}
        Lemp.boot <- function(i) Emp.boot(i, i)/(i/n)
        boot.res <- rbind(boot.res, Vectorize(Lemp.boot)(1:(n-1)))
      }
      lm <- apply(boot.res, 2, mean)
      lq <- apply(boot.res[,1:(n-1)], 2, stats::quantile, .025)
      uq <- apply(boot.res[,1:(n-1)], 2, stats::quantile, .975)
      if (!control$smooth){
        wm <- which.min(lm[((n-1)*(1-control$cutoff)):(n-1)]) + (n-1)*(1-control$cutoff)
        m  <- min(lm[((n-1)*(1-control$cutoff)):(n-1)])
      }
      if (control$smooth) {
        con <- 1:(n-1)
        emp.lamb.ss <- stats::predict(stats::loess(lm ~ con, span = control$smoothing.span))
        wm <- which.min(emp.lamb.ss[((n-1)*(1-control$cutoff)):(n-1)]) + (n-1)*(1-control$cutoff)
        m <- min(emp.lamb.ss[((n-1)*(1-control$cutoff)):(n-1)])
      }
      if (plot){
        plot(boot.res[1,], ylim=c(0,1), type="l", ylab="lower tail dependence")
        for (k in 2:control$bootstrap.num) lines(boot.res[k,])
        lines(lm, col=2, lwd=2)
        lines(lq, col=2)
        lines(uq, col=2)
        if (control$smooth) lines(emp.lamb.ss, col=4, lwd=2)
        abline(h=m, v=wm)
      }
    }
    return(c(m, wm))
  }
}

