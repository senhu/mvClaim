#' @title Density of bivariate gamma distribution given parameters values
#' @description Density of bivariate gamma distribution given parameters values
#'
#' @param y vector of quantiles
#' @param alpha vector of alpha values in bivariate gamma distribution
#' @param beta value of beta in bivariate gamma distribution
#' @param log logical; if TRUE, probabilities p are given as log(p). Default is FALSE.
#'
#' @export
#'
#' @example
#' dbivgamma(c(1,1), alpha=c(.3,.4,.5), beta=1, log = TRUE)
#'
#'
#' c1 <- seq(0, 15, by=0.1)[-1]
#' c2 <- seq(0, 15, by=0.1)[-1]
#' restable <- matrix(0, nrow = length(c1), ncol=length(c2))
#' for (i in c(1:length(c1))){
#'  for (j in c(1:length(c2))){
#'    restable[i,j] <- dbivgamma(c(c1[i], c2[j]),
#'                               alpha = c(0.5,0.5,0.5),
#'                               beta = 1,log=TRUE)
#'  }
#' }
#' rownames(restable) <- NULL
#' colnames(restable) <- NULL
#' require(lattice)
#' new.palette=colorRampPalette(c("white","yellow","red"),space="rgb")
#' levelplot(restable, col.regions=new.palette(20), xlab="Y_1", ylab="Y_2")

dbivgamma <- function(y, alpha, beta, log=FALSE){
  y1 <- y[1]
  y2 <- y[2]
  alpha1 <- alpha[1]
  alpha2 <- alpha[2]
  alpha3 <- alpha[3]
  beta <- beta
  integrand <- function(x3){
    dgamma(y1-x3,shape=alpha1,rate=beta, log=FALSE)*
      dgamma(y2-x3,shape=alpha2,rate=beta, log=FALSE)*
      dgamma(x3,shape=alpha3,rate=beta, log=FALSE)
  }
  log.integrand <- function(x3){
    dgamma(y1-x3,shape=alpha1,rate=beta, log=TRUE)+
      dgamma(y2-x3,shape=alpha2,rate=beta, log=TRUE)+
      dgamma(x3,shape=alpha3,rate=beta, log=TRUE)
  }
  int <- distrEx::distrExIntegrate(integrand, lower=0, upper=min(y1,y2),
                                   rel.tol = 1e-6, order=1e+3)

  if (!isTRUE(log) && int==Inf){ return(list("value"=int, "message"="int is infinity, should use log==TRUE"))}
  if (!isTRUE(log) && int!=Inf){ return(list("value"=int, "message"="OK")) }

  if (isTRUE(log) && int!=Inf && log(int)!=-Inf ){return(list("value"=log(int), "message"="OK"))}
  #if (){return(log(int)); warning("int is 0, should use log==FALSE") }
  if ((isTRUE(log) && int==Inf) || (isTRUE(log) && int!=Inf && log(int)==-Inf)){
    #log.int <- myLogIntegrate(log.integrand, lower=.Machine$double.eps^{1/3}, upper=(min(y1, y2)-.Machine$double.eps^{1/3}), tol=1, verbose=FALSE)
    #log.int <- log.simpson(log.integrand, a=0, b=min(y1,y2), n=100000)
    log.int <- LogSimpson(log.integrand,
                           lower=.Machine$double.eps^(1/3),
                           upper=min(y1,y2)-.Machine$double.eps^(1/3),
                           n=1e+5)
    if (abs(log.int)!=Inf) {return(list("value"=log.int, "message"="using log-integral approximation"))} else {
      stop("error in dbivgamma: logIntegrate function")
    }}
}


LogSimpson <- function(logfun, lower, upper, n=100) {
  # numerical integral using Simpson's rule
  # assume lower < upper and n is an even positive integer
  h <- (upper-lower)/n
  x <- seq(lower, upper, by=h)
  if (n == 2) {
    s <- matrixStats::logSumExp(c(logfun(x[1]), logfun(x[2]), logfun(x[2]), logfun(x[2]), logfun(x[2]), logfun(x[3])))
  } else {
    s <- matrixStats::logSumExp(c(logfun(x[1]),
                                  logfun(x[n+1]),
                                  logfun(x[seq(2,n,by=2)]),
                                  logfun(x[seq(2,n,by=2)]),
                                  logfun(x[seq(3,n-1, by=2)]),
                                  logfun(x[seq(3,n-1, by=2)]),
                                  logfun(x[seq(3,n-1, by=2)]),
                                  logfun(x[seq(3,n-1, by=2)])))
  }
  s <- log(h)-log(3) + s
  return(s)
}