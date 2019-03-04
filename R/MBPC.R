#' Model-based bivariate Poisson clustering
#'
#' @param G number of mixtures
#' @param data
#' @param maxit maximum iteration. Default is 300.
#' @param tol tolerance. Default is 1e-5.
#' @param Aitken if TRUE, using Aitken acceleration estimate criterion. Defaul is TRUE.
#' @param verbose if TRUE, print the iteration results. Default is TRUE.


MBPC <- function(G,
                 data,
                 start,
                 maxit=300,
                 tol=1e-5,
                 Aitken=FALSE,
                 verbose=TRUE){

  options(warn=-1) 
  y1<-data[,1]
  y2<-data[,2]
  n<-length(y1)
  p<-dim(as.data.frame(data))[2]
  if (p!=2){stop("The data set is not bivariate.")}
  
  switch(start, 
         random.start={
           index <- sample(c(1:G), n, replace=TRUE)
           p.z <- table(index)/n
           lambda1 <- rep(0,G)
           lambda2 <- rep(0,G)
           lambda3 <- rep(0,G)
           for (g in seq_len(G)){
             temp.dat <- data[index==g,]
             temp.v1 <- MASS::fitdistr(temp.dat[,1], "poisson")
             temp.v2 <- MASS::fitdistr(temp.dat[,2], "poisson")
             lambda1.temp <- temp.v1$estimate[1]
             lambda2.temp <- temp.v2$estimate[1]
             lambda3[g] <- min(lambda1.temp, lambda2.temp)/3
             lambda1[g] <- lambda1.temp - lambda3[g]
             lambda2[g] <- lambda2.temp - lambda3[g]
           }
         },
         kmeans.start={
           
         },
         hierarchical.start={
           
         },
         given.start={
           lambda1<-c(0.01, 1)
           lambda2<-c(0.01, 1)
           lambda3<-c(0.01, 1)
         })

  zero <- ( y1==0 )|( y2==0 )
  j<-1
  loglike.diff <- 1000
  loglike.current <- -sqrt(.Machine$double.xmax)
  loglike<-rep(0,maxit)
  loglike.store <- rep(-sqrt(.Machine$double.xmax), 3)
  lambda1.track <- NULL
  lambda2.track <- NULL
  lambda3.track <- NULL
  while ( (loglike.diff > tol) && (j <= maxit) ) {
    #--------
    #E step  
    #--------
    den <- matrix(0, nrow=n, ncol=G)
    pzy <- NULL
    s<-matrix(0,nrow=n, ncol=G)
    for (i in 1:n) {
      for (g in seq_len(G)) {
        den[i,g] <- dbivpois(c(y1[i],y2[i]),
                             lambda=c(lambda1[g],lambda2[g],lambda3[g]),
                             log=TRUE)
        if (zero[i]) { s[i, g] <- 0 }
        else {
          ll.alt <-dbivpois(c(y1[i]-1,y2[i]-1),
                            lambda=c(lambda1[g],lambda2[g],lambda3[g]),
                            log=TRUE)
          s[i,g] <- exp(log(lambda3[g]) + ll.alt - den[i,g])
          if (is.nan(s[i,g]) || is.na(s[i,g])){s[i, g]<-0}
        }
      }
    }
    newden <- sweep(den, 2, log(p.z), FUN="+")
    denom  <- matrixStats::rowLogSumExps(newden)
    pzy    <- exp(sweep(newden, 1, denom, FUN="-"))
    #if (table(round(rowSums(pzy),2))!=n){stop("Wrong E-step estimation: pzy")}

    loglike.new <- sum(denom)
    loglike[j] <- loglike.new
    loglike.store <- c(loglike.store[-1], loglike.new)
    if (Aitken){
      loglike.diff <- abs( MoEClust::MoE_aitken(loglike.store)$linf - loglike.current )
    } else {
      loglike.diff <- abs( loglike.new - loglike.current ) / (1+abs(loglike.new))
    }
    loglike.current<-loglike.new
    
    #--------
    # M step  
    #--------
    p.z <- colSums(pzy)/n
    for (g in seq_len(G)){
      x1<-y1-s[,g]
      x2<-y2-s[,g]
      x3<-s[,g]
      lambda1[g]<- sum(pzy[,g]*x1)/sum(pzy[,g])
      lambda2[g]<- sum(pzy[,g]*x2)/sum(pzy[,g])
      lambda3[g]<- sum(pzy[,g]*x3)/sum(pzy[,g])
    }
    
    if (verbose){
      cat(c("iter:", j, "; current loglike:", round(loglike.new,6), "; loglike change:", round(loglike.diff,6), "\n"))
    }
    
    lambda1.track <- rbind(lambda1.track, lambda1)
    lambda2.track <- rbind(lambda2.track, lambda2)
    lambda3.track <- rbind(lambda3.track, lambda3)
    
    j<-j+1
  }

  noparams<- G*4 - 1
  AIC<- -2*loglike[j-1] + noparams * 2
  BIC<- -2*loglike[j-1] + noparams * log(n)

  classification <- apply(pzy, 1, which.max)
  
  result<-list(lambda1=lambda1,
               lambda2=lambda2,
               lambda3=lambda3,
               pro = p.z,
               pzy = pzy,
               loglikelihood=loglike[j-1],
               loglike.seq=loglike[1:(j-1)],
               classification=classification,
               parameters.number=noparams,
               AIC=AIC,
               BIC=BIC,
               iterations=(j-1))
  options(warn=0)
  class(result)<-c('MBPC')
  return(result)
}
