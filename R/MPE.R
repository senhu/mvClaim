#' Multivariate Poisson distribution estimation
#'

bvpe <- function(data, method="EM", verbose=TRUE, tol=1e-6){
  n <- dim(data)[1]
  p <- dim(data)[2]
  if (p!=2){stop("not bivariate")}
  y1 <- data[,1]
  y2 <- data[,2]
  zero <- ( y1==0 )|( y2==0 )
  maxit=1e+4
  loglikelihood <- rep(0, maxit)
  lambda1 <- 1
  lambda2 <- 4
  lambda3 <- 2
  iter = 1
  loglike.diff <- 1000
  loglike.current <- -.Machine$double.eps
  den <- rep(0, n)
  s<-rep(0, n)
  while ( (loglike.diff > tol) && (iter <= maxit) ){

    for (i in 1:n) {
      den[i] <- dbivpois(c(y1[i],y2[i]),
                         lambda=c(lambda1,lambda2,lambda3),
                         log=TRUE)
      if (zero[i]) { s[i] <- 0 }
      else {
        ll.alt <-dbivpois(c(y1[i]-1,y2[i]-1),
                          lambda=c(lambda1,lambda2,lambda3),
                          log=TRUE)
        s[i] <- exp(log(lambda3) + ll.alt - den[i])
        if (is.nan(s[i]) || is.na(s[i])){s[i, g]<-0}
      }
    }

    lambda1 <- sum(y1-s)/n
    lambda2 <- sum(y2-s)/n
    lambda3 <- sum(s)/n

    loglike.new <- sum(den)
    loglikelihood[iter] <- loglike.new
    loglike.diff <- abs( loglike.new - loglike.current ) / (1+abs(loglike.new))
    if (verbose){
      cat(c("iter:", iter, "; current loglike:", round(loglike.new,6), "; loglike change:", round(loglike.diff,6), "\n"))
    }
    iter <- iter + 1
    loglike.current <- loglike.new
  }
  return(parameter=c(lambda1, lambda2, lambda3))
}
