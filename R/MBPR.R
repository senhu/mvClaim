#' Mixture of bivariate Poisson regression (mixture of experts)
#'
#' @param K number of mixtures
#' @param l1 formula for lambda1
#' @param l2 formula for lambda2
#' @param l3 formula for lambda3. If NULL, only double Poisson regression. Default is NULL.
#' @param lp formula for mixing proportion. If NULL, mixing proportion is not regressed against covariates. Default is NULL.
#' @param w weights in GLM (if applicable) such as "exposure". Default is NULL.
#' @param data
#' @param maxit maximum iteration. Default is 300.
#' @param tol tolerance. Default is 1e-5.
#' @param Aitken if TRUE, using Aitken acceleration estimate criterion. Defaul is TRUE.
#' @param verbose if TRUE, print the iteration results. Default is TRUE.


MBPR<- function(G,
                data,
                l1,
                l2,
                l3,
                lp=NULL,
                w=NULL,
                maxit=300,
                tol=1e-6,
                Aitken=FALSE,
                verbose=TRUE){

  options(warn=-1)
  templist <- list(G=G,
  						l1=l1,
  						l2=l2,
  						l3=l3,
  						lp=lp,
  						w =w,
  						data=substitute(data),
  						maxit=maxit,
  						tol=tol,
  						Aitken=Aitken,
  						verbose=verbose)
  tempcall<-as.call( c(expression(MBPR), templist))
  rm(templist)

  if (is.null(l3)) {stop("l3 is NULL")}
  
  namey1 <- as.character(l1[2])
  namey2 <- as.character(l2[2])
  y1     <- data[,names(data)==namey1]
  y2     <- data[,names(data)==namey2]
  n      <- length(y1)
  #p      <- dim(as.data.frame(data))[2]
  if (!is.null(w)) {weight.var <- data[, names(data)==w]} else {weight.var <- rep(1, n)}
  data$weight.var <- weight.var
  data1 <- data; data1 <- data1[,names(data1)!=namey1]; data1 <- data1[ , names(data1)!=namey2]
  data2<-data1; data3<-data1
  datap<-data1
  
  zero <- ( y1==0 )|( y2==0 )

  if (!is.null(l1)){
    if (as.character(l1[3])=="."){l1.new <- formula( paste( "x1", paste( names(data1),'',collapse='+',sep='' ), sep='~'))}
    else {l1.new <- formula(paste("x1", as.character(l1)[3], sep="~"))}
  }
  if (!is.null(l2)){
    if (as.character(l2[3])=="."){l2.new <- formula( paste( "x2", paste( names(data1),'',collapse='+',sep='' ), sep='~'))}
    else {l2.new <- formula(paste("x2", as.character(l2)[3], sep="~"))}
  }
  if (!is.null(l3)){
    if (as.character(l3[2])=="."){l3.new <- formula( paste( "x3", paste( names(data1), "", collapse = "+", sep = ""), sep="~"))}
    else {l3.new <- formula(paste("x3", as.character(l3)[2], sep="~"))}
  }
  if (!is.null(lp)){
    if (as.character(lp[2])=="."){lp.new <- formula( paste( "pzy", paste( names(data1), "", collapse = "+", sep = ""), sep="~"))}
    else {lp.new <- formula(paste("pzy", as.character(lp)[2], sep="~"))}
  }

  # starting values
  if (is.null(lp)){
  	tm <- runif(G); p.z <- tm/sum(tm)
  } else {
    tm <- runif(G); p.z.row <- tm/sum(tm); p.z <- matrix(rep(p.z.row,each=n), nrow = n, byrow=FALSE)
  }
  lambda1.temp    <- glm(l1, family=poisson, data=data)$fitted
  lambda1.current <- lambda1.temp * t(rmultinom(n, size=100, prob = rep(1/G, G))/100)
  lambda2.temp    <- glm(l2, family=poisson, data=data)$fitted
  lambda2.current <- lambda2.temp * t(rmultinom(n, size=100, prob = rep(1/G, G))/100)
  lambda3.temp    <- rep( max(0.1, cov(y1,y2,use='complete.obs')), n)
  lambda3.current <- lambda3.temp * t(rmultinom(n, size=100, prob = rep(1/G, G))/100)
  
  
  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.xmax)
  j               <- 1
  loglike         <- rep(0,maxit)
  lambda1.trace   <- NULL
  lambda2.trace   <- NULL
  lambda3.trace   <- NULL
  coef1.trace     <- NULL
  coef2.trace     <- NULL
  coef3.trace     <- NULL
  coefp.trace     <- NULL
  loglike.store   <- rep(-.Machine$double.xmax, 3)
    
  while ( (loglike.diff > tol) && (j <= maxit) ) {
    #--------
    # E step  
    #--------
    s   <- matrix(0,nrow=n, ncol=G)
  	den <- matrix(0, nrow=n, ncol=G)
    for (i in 1:n) {
   		for (g in 1:G) {
   		  den[i,g] <- dbivpois(c(y1[i],y2[i])  ,
                              lambda=c(lambda1.current[i, g],
                                       lambda2.current[i, g],
                                       lambda3.current[i, g]),
                              log=TRUE)
          if (zero[i]) { s[i, g] <- 0 }
          else {
            den.alt <- dbivpois(c(y1[i]-1,y2[i]-1),
                                lambda=c(lambda1.current[i, g],
                                         lambda2.current[i, g],
                                         lambda3.current[i, g]),
                                log=TRUE)
            s[i,g] <- exp( log(lambda3.current[i,g]) + den.alt - den[i,g] )
            if (is.nan(s[i,g]) || is.nan(s[i,g])){s[i,g]<-0; warning("warning1: E-step: s")}
          }
        }
      }
      if (is.null(lp)){newden <- sweep(den, 2, log(p.z), FUN="+")} else{
        newden <- den + log(p.z)
      }
      denom  <- matrixStats::rowLogSumExps(newden)
      z.new    <- exp(sweep(newden, 1, denom, FUN="-"))
      #if (table(round(rowSums(z.new),2))!=n){stop("Wrong E-step estimation: pzy")}
      
      loglike.new    <- sum(denom)
      loglike[j]     <- loglike.new
      loglike.store  <- c(loglike.store[-1], loglike.new)
      if (Aitken){
        loglike.diff <- abs( MoEClust::MoE_aitken(loglike.store)$linf - loglike.current )
      } else {
        loglike.diff <- abs( loglike.new - loglike.current ) / (1+abs(loglike.new))
      }
      loglike.current<-loglike.new
      
      #--------
      # M step
      #--------
      if (is.null(lp)){ p.z <- colSums(z.new)/n } else {
        datap$pzy <- z.new
        mp <- nnet::multinom(formula = lp.new, data = datap, trace=FALSE) # reltol
        coef.p <- coef(mp)
        p.z <- stats::predict(mp, type="probs")
      }

      coef1 <- coef2 <- coef3 <- NULL
      lambda1.new <- lambda2.new <- lambda3.new <- matrix(0,nrow=n, ncol=G)
      for (g in c(1:G)){
        x1 <- abs(y1-s[,g]); if (min(x1)<0) warning("min of x1 negative")
        x2 <- abs(y2-s[,g]); if (min(x2)<0) warning("min of x2 negative")
        x3 <- abs(s[,g])
        
        data1$x1        <- x1
        m1              <- glm(l1.new,family=poisson(link="log"),data=data1,weights=(z.new[,g]))
        coef1           <- cbind(coef1, m1$coef)
        lambda1.new[,g] <- m1$fitted
        
        data1$x2        <- x2
        m2              <- glm(l2.new,family=poisson(link="log"),data=data1,weights=(z.new[,g]))
        coef2           <- cbind(coef2, m2$coef)
        lambda2.new[,g] <- m2$fitted
        
        data1$x3        <- x3
        m3              <- glm(l3.new,family=poisson(link="log"),data=data1,weights=(z.new[,g]))
        coef3           <- cbind(coef3, m3$coef)
        lambda3.new[,g] <- m3$fitted
      }
      if (verbose){
        cat(c("iter:", j, "; current loglike:", round(loglike.new,6), "; loglike change:", round(loglike.diff,6), "\n"))
      }

      lambda1.trace   <- append(lambda1.trace, list(lambda1.new))
      lambda2.trace   <- append(lambda2.trace, list(lambda2.new))
      lambda3.trace   <- append(lambda3.trace, list(lambda3.new))
      coef1.trace     <- append(coef1.trace, list(coef1))
      coef2.trace     <- append(coef2.trace, list(coef2))
      coef3.trace     <- append(coef3.trace, list(coef3))
      lambda1.current <- lambda1.new
      lambda2.current <- lambda2.new
      lambda3.current <- lambda3.new
      
      j               <- j+1
    }

    if (is.null(lp)){
      noparams<- G*(m1$rank + m2$rank + m3$rank) + (G-1)
    } else {
      noparams<- G*(m1$rank + m2$rank + m3$rank) + mp$rank
    }
    AIC<- -2*loglike[j-1] + noparams * 2
    BIC<- -2*loglike[j-1] + noparams * log(n)
    
    if (is.null(lp)){
      y1.fitted <- (lambda1.current + lambda3.current) %*% p.z
      y2.fitted <- (lambda2.current + lambda3.current) %*% p.z
    } else {
      y1.fitted <- rowSums(p.z * (lambda1.current + lambda3.current))
      y2.fitted <- rowSums(p.z * (lambda2.current + lambda3.current))
    }
    

    classification <- apply(z.new, 1, which.max)

    #  Calculation of output
    result<-list(coefficients  = c(list(coef1), list(coef2), list(coef3)),
                 lambda1       = lambda1.current,
                 lambda2       = lambda2.current,
                 lambda3       = lambda3.current,
                 pro           = p.z,
                 z             = z.new,
                 fitted.values = cbind(y1.fitted,
                                       y2.fitted),
                 residuals     = cbind(y1=(y1-y1.fitted),
                                       y2=(y2-y2.fitted)),
                 loglikelihood = loglike[j-1],
                 loglike.seq   = loglike[1:(j-1)],
                 class         = classification,
                 parameters.no = noparams,
                 AIC           = AIC,
                 BIC           = BIC,
                 iterations    = (j-1),
                 call          = tempcall)
                 
  options(warn=0)
  class(result)<-c('MBPR', 'glm')
  return(result)
}
