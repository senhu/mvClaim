rbivpois <- function(n, lambda){
  lambda1 <- lambda[1]
  lambda2 <- lambda[2]
  lambda3 <- lambda[3]
  x1 <- rpois(n, lambda = lambda1)
  x2 <- rpois(n, lambda = lambda1)
  x3 <- rpois(n, lambda = lambda1)
  y1 <- x1+x3
  y2 <- x2+x3  
  return(cbind(y1,y2))
}

ha <- rbivpois(800, lambda=c(0.93307588, 3.87819783, 0.04253387))
plot(jitter(ha))
