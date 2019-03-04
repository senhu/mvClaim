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


#haha1 <- rbivgamma(2000, c(1,1,1),beta=1)
#plot(haha1)
# haha2 <- rbivgamma(2000, c(0.5406931,0.5399668,0.2699834),beta=0.002210341)
# plot(rbind(haha1, haha2))
# plot(data)
