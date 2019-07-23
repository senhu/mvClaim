#' Density of Bivariate Gamma Distribution
#'
#' Density of bivariate gamma distribution given parameters values
#'
#' @param y vector of quantiles
#' @param alpha vector of alpha values in bivariate gamma distribution
#' @param beta value of beta in bivariate gamma distribution
#' @param log logical; if TRUE, probabilities p are given as log(p). Default is FALSE.
#'
#' @examples
#' dbivgamma(c(1,1), alpha=c(0.3, 0.4, 0.5), beta=1, log = TRUE)
#'
#'
#' @importFrom stats Gamma dgamma
#' @export

dbivgamma <- function(y, alpha, beta, log=FALSE){
  y1 <- y[1]
  y2 <- y[2]
  alpha1 <- alpha[1]
  alpha2 <- alpha[2]
  alpha3 <- alpha[3]
  beta <- beta
  integrand <- function(x3){
    stats::dgamma(y1-x3,shape=alpha1,rate=beta, log=FALSE)*
      stats::dgamma(y2-x3,shape=alpha2,rate=beta, log=FALSE)*
      stats::dgamma(x3,shape=alpha3,rate=beta, log=FALSE)
  }
  log.integrand <- function(x3){
    stats::dgamma(y1-x3,shape=alpha1,rate=beta, log=TRUE)+
      stats::dgamma(y2-x3,shape=alpha2,rate=beta, log=TRUE)+
      stats::dgamma(x3,shape=alpha3,rate=beta, log=TRUE)
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
