rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
rep.row <- function(x,n){
  matrix(rep(x,each=n), nrow = n, byrow=FALSE)
}
# rep.row(c(1,2,3),5)
# rep.col(c(1,2,3),5)
