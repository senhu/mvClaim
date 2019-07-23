
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
rep.row <- function(x,n){
  matrix(rep(x,each=n), nrow = n, byrow=FALSE)
}
# rep.row(c(1,2,3),5)
# rep.col(c(1,2,3),5)


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
