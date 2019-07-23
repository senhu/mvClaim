#' Random generation for the bivariate gamma distribution
#'
#' This function allows you to generate random bivariate gamma variables given the distribution parameters
#' @param n number of observations.
#' @param alpha the shape parameters of the bivariate gamma distribution. Must be a vector of length 3 represent alpha 1, alpha 2 and alpha 3 respectively.
#' @param beta the rate parameter of the bivariate gamma distribution.
#'
#' @return random deviates of length n
#'
#' @examples
#' plot(rbivgamma(2000, alpha=c(1,1,1), beta=1))
#' plot(rbivgamma(2000, c(0.5, 0.5, 0.2),beta=0.001))

rbivgamma <- function(n, alpha, beta){
  alpha1 <- alpha[1]
  alpha2 <- alpha[2]
  alpha3 <- alpha[3]
  g1 <- rgamma(n, alpha1, beta)
  g2 <- rgamma(n, alpha2, beta)
  g3 <- rgamma(n, alpha3, beta)
  y1 <- g1+g3; y2 <- g2+g3
  return(cbind(y1,y2))
}



