#' Density of bivariate Poisson distribution
#'
#' @param y a vector
#' @param lambda vector of lambda values for lambda1, lambda2, lambda3.
#' @param method "analytical" or "recursive"
#' @param log whether the calculated density is on log scale
#' @note This function is borrowed from "bivpois" R package which is no longer available on CRAN
#' @example
#' dbivpois(y=c(2, 2), lambda = c(1,1,1), method = "recursive")
#' dbivpois(y=c(2, 2), lambda = c(1,1,1), method = "analytical")
#' dmvpois(y=c(2,2), lambda=c(1,1,1), method="recursive", log = F)
#'
#' library(lattice)
#' new.palette=colorRampPalette(c("white","yellow","red"),space="rgb")
#' tab1 <- dbivpois(y=c(15, 15), lambda =c(2, 3, 1), method="recursive")
#' rownames(tab1) <- c(0:15)
#' colnames(tab1) <- c(0:15)
#' levelplot(tab1, col.regions=new.palette(20), xlab="Y_1", ylab="Y_2")

dbivpois <- function(y, lambda = c(1, 1, 1), log=FALSE, method="analytical") {
  if (!is.vector(y)){
    stop("y is not a vector")
  }
  if (is.vector(y) && length(y) != 2) {
    stop("Bivariate Poisson is only defined for vectors of length 2")}
  y1<-y[1]
  y2<-y[2]
  lambda1 <- lambda[1]
  lambda2 <- lambda[2]
  lambda3 <- lambda[3]
  switch(method,
         recursive={
           if(sum(y == 0)==2){
             prob <- exp(-sum(lambda))
           }
           else {
             prob <- matrix(NA, nrow = y1+1, ncol = y2+1, byrow = T)
             prob[1, 1] <- exp( - lambda1 - lambda2 - lambda3)
             for(i in 2:(y1 + 1)) {
               prob[i, 1] <- (prob[i - 1, 1] * lambda1)/(i - 1)
             }
             for(j in 2:(y2 + 1)) {
               prob[1, j] <- (prob[1, j - 1] * lambda2)/(j - 1)
             }
             for(j in 2:(y2 + 1)) {
               for(i in 2:(y1 + 1)) {
                 prob[i, j] <- (lambda1 * prob[i-1, j]+lambda3 * prob[i-1, j-1])/(i - 1)
               }
             }
           }
           result <- prob
         },
         analytical={
             xymin<-min( y1,y2 )
             lambdaratio<-lambda3/(lambda1*lambda2)

             i <- 0:xymin
             sums<- -lgamma(y1-i+1)-
                    lgamma(i+1)-
                    lgamma(y2-i+1)+
                    i*log(lambdaratio)

             maxsums <- max(sums)
             sums<- sums - maxsums
             logsummation<- log( sum(exp(sums)) ) + maxsums
             logbp <- -sum(lambda) + y1 * log( lambda[1] ) + y2 * log( lambda[2] ) + logsummation

           if (log) { result<-    logbp }
           else     { result<-exp(logbp)  }
         },
         stop('Unknown method. Please use "recursive", or "analytical"'))
  return(result)
}



