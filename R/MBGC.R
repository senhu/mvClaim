#' Mixture of bivariate gamma distributions clustering
#'
#' Estimation using EM algorithm for mixture of bivariate gamma distributions



MBGC <- function(data,
                 G,
                 maxit = 300,
                 tol = 1e-6,
                 start = "random.start",
                 start.value=NULL,
                 Aitken = FALSE,
                 verbose= TRUE){

  #options(warn=-1)
  templist <- list(data   = substitute(data),
                   G      = G,
                   maxit  = maxit,
                   tol    = tol,
                   start  = start,
                   start.value = start.value,
                   Aitken = Aitken,
                   verbose= verbose)
  tempcall <- as.call( c(expression(MBGC), templist) ); rm(templist)

  n<-dim(data)[1]
  p<-dim(data)[2]
  if (p!=2){stop("The data set is not bivariate.")}

  if (start=="random.start"){
    index          <- sample(c(1:G), n, replace=TRUE, prob=c(1/3,rep((2/3)/(G-1), G-1)))
    p.z.current    <- table(index)/n
    alpha1.current <- rep(0,G)
    alpha2.current <- rep(0,G)
    alpha3.current <- rep(0,G)
    beta.current   <- rep(0,G)
    for (g in seq_len(G)){
      temp.dat <- data[index==g,]
      temp.v1 <- MASS::fitdistr(temp.dat[,1]/mean(temp.dat[,1]), "gamma")
      temp.v2 <- MASS::fitdistr(temp.dat[,2]/mean(temp.dat[,2]), "gamma")
      beta.current[g] <- mean(c(temp.v1$estimate[2]/mean(temp.dat[,1]),
                                temp.v2$estimate[2]/mean(temp.dat[,2])))
      alpha1.temp <- temp.v1$estimate[1]
      alpha2.temp <- temp.v2$estimate[1]
      alpha3.current[g] <- min(alpha1.temp, alpha2.temp)/3
      alpha1.current[g] <- alpha1.temp - alpha3.current[g]
      alpha2.current[g] <- alpha2.temp - alpha3.current[g]
    }
  }
  if (start=="mclust"){

    initial.mclust <- mclust::Mclust(G=G, data = cbind(y1,y2))
    p.z.current        <- initial.mclust$parameters$pro
    index          <- initial.mclust$classification
    alpha1.current <- rep(0,G)
    alpha2.current <- rep(0,G)
    alpha3.current <- rep(0,G)
    beta.current   <- rep(0,G)
    for (g in seq_len(G)){
      temp.dat <- data[index==g,]
      temp.v1 <- MASS::fitdistr(temp.dat[,1]/mean(temp.dat[,1]), "gamma")
      temp.v2 <- MASS::fitdistr(temp.dat[,2]/mean(temp.dat[,2]), "gamma")
      beta.current[g] <- mean(c(temp.v1$estimate[2]/mean(temp.dat[,1]),
                                temp.v2$estimate[2]/mean(temp.dat[,2])))
      alpha1.temp <- temp.v1$estimate[1]
      alpha2.temp <- temp.v2$estimate[1]
      alpha3.current[g] <- min(alpha1.temp, alpha2.temp)/3
      alpha1.current[g] <- alpha1.temp - alpha3.current[g]
      alpha2.current[g] <- alpha2.temp - alpha3.current[g]
    }
  }
  if (start=="given"){
    alpha1.current=start.value[1,]
    alpha2.current=start.value[2,]
    alpha3.current=start.value[3,]
    beta.current=start.value[4,]
    p.z.current=start.value[5,]
  }

  j               <- 1
  loglike.diff    <- 1000
  loglike         <- rep(0,maxit)
  loglike.store   <- rep(-sqrt(.Machine$double.xmax), 3)
  alpha1.trace    <- alpha1.current
  alpha2.trace    <- alpha2.current
  alpha3.trace    <- alpha3.current
  beta.trace      <- beta.current
  #z.trace       <- NULL
  loglike.current <- -sqrt(.Machine$double.xmax)

  while ( (loglike.diff >= tol) && (j <= maxit) ) {
    #-------
    #E step
    #-------
    Expected.x3    <- matrix(0,nrow=n, ncol=G)
    Expected.x1    <- matrix(0,nrow=n, ncol=G)
    Expected.x2    <- matrix(0,nrow=n, ncol=G)
    Expected.logx3 <- matrix(0,nrow=n, ncol=G)
    Expected.logx1 <- matrix(0,nrow=n, ncol=G)
    Expected.logx2 <- matrix(0,nrow=n, ncol=G)
    for (g in seq_len(G)){
      Estep.res <- MBGC.Estep.function(data,
                                       alpha=c(alpha1.current[g], alpha2.current[g], alpha3.current[g]),
                                       beta=beta.current[g])
      Expected.x3[,g]    <- Estep.res[[1]]
      Expected.x1[,g]    <- (data[,1]-Estep.res[[1]])
      Expected.x2[,g]    <- (data[,2]-Estep.res[[1]])
      Expected.logx3[,g] <- Estep.res[[2]]
      Expected.logx1[,g] <- Estep.res[[3]]
      Expected.logx2[,g] <- Estep.res[[4]]
    }

    den <- matrix(0, nrow=n, ncol=G)
    for (i in seq_len(n)) {
      for (g in seq_len(G)) {
        den[i,g] <- dbivgamma(c(data[i,1],data[i,2]) ,
                              alpha=c(alpha1.current[g],alpha2.current[g],alpha3.current[g]),
                              beta=beta.current[g],
                              log=TRUE)$value
        if (den[i,g]==-Inf){stop("warings1")}
      }
    }
    newden <- sweep(den, 2, log(p.z.current), FUN="+")
    denom  <- matrixStats::rowLogSumExps(newden)
    z.new    <- exp(sweep(newden, 1, denom, FUN="-"))
    # z.trace <- append(z.trace, list(z.new))
    if (table(round(rowSums(z.new),0))!=n){stop("Wrong E-step estimation: z")}

    #-------
    #M step
    #-------
    alpha1.new <- rep(0, G)
    alpha2.new <- rep(0, G)
    alpha3.new <- rep(0, G)
    beta.new   <- rep(0, G)

    p.z.new    <- colSums(z.new)/n

    for (g in seq_len(G)){
      beta.new[g] <- ((alpha1.current[g]+alpha2.current[g]+alpha3.current[g])*sum(z.new[,g])) / (z.new[,g]%*%(Expected.x1[,g]+Expected.x2[,g]+Expected.x3[,g]))
      alpha1.rootfun <- function(alpha1.var,wh){
        log(beta.new[wh])*sum(z.new[,wh]) - digamma(alpha1.var)*sum(z.new[,wh]) + z.new[,wh]%*%Expected.logx1[,g]
      }
      alpha1.new[g] <- uniroot(alpha1.rootfun,wh=g,
                               lower=.Machine$double.eps, #sqrt(.Machine$double.eps),
                               upper=100000,
                               tol = sqrt(.Machine$double.xmin))$root
      alpha2.rootfun <- function(alpha2.var,wh){
        log(beta.new[wh])*sum(z.new[,wh]) - digamma(alpha2.var)*sum(z.new[,wh]) + z.new[,wh]%*%Expected.logx2[,g]
      }
      alpha2.new[g] <- uniroot(alpha2.rootfun,wh=g,
                               lower=.Machine$double.eps,
                               upper=100000,
                               tol = sqrt(.Machine$double.xmin))$root
      alpha3.rootfun <- function(alpha3.var,wh){
        log(beta.new[wh])*sum(z.new[,wh]) - digamma(alpha3.var)*sum(z.new[,wh]) + z.new[,wh]%*%Expected.logx3[,g]
      }
      alpha3.new[g] <- uniroot(alpha3.rootfun,wh=g,
                               lower=.Machine$double.eps,
                               upper=100000,
                               tol=sqrt(.Machine$double.xmin))$root
    }

    if (verbose){
      cat(c("alpha1:", round(alpha1.new,4), "\n"))
      cat(c("alpha2:", round(alpha2.new,4), "\n"))
      cat(c("alpha3:", round(alpha3.new,4), "\n"))
      cat(c("beta:", round(beta.new,4), "\n"))
      cat(c("p.z:", round(p.z.new,4), "\n"))
    }

    loglike.new <- sum(denom)
    loglike[j]  <- loglike.new
    loglike.store <- c(loglike.store[-1], loglike.new)
    if (Aitken){
      loglike.diff <- abs(MoEClust::MoE_aitken(loglike.store)$linf - loglike.current)
    } else {
      loglike.diff <- abs(loglike.new - loglike.current) / (1+abs(loglike.new))
    }

    if (verbose){
      cat(c("iter:", j, "; current loglike:", round(loglike.new,4), "; loglike change:", round(loglike.diff,7), "\n"))
    }

    alpha1.trace   <- rbind(alpha1.trace, alpha1.new)
    alpha2.trace   <- rbind(alpha2.trace, alpha2.new)
    alpha3.trace   <- rbind(alpha3.trace, alpha3.new)
    beta.trace     <- rbind(beta.trace, beta.new)
    alpha1.current <- alpha1.new
    alpha2.current <- alpha2.new
    alpha3.current <- alpha3.new
    beta.current   <- beta.new
    p.z.current    <- p.z.new
    loglike.current<- loglike.new
    j              <- j+1
  }

  if (!all(loglike[1:(j-1)] == cummax(loglike[1:(j-1)]))){
    warning("Loglikelihood is not monotonically increasing!")
  }
  classification <- apply(z.new, 1, which.max)
  noparams       <- 5*G - 1
  AIC            <- -2*loglike.new + noparams * 2
  BIC            <- -2*loglike.new + noparams * log(n)

  result<-list(alpha1        = alpha1.current,
               alpha2        = alpha2.current,
               alpha3        = alpha3.current,
               beta          = beta.current,
               pro           = p.z.current,
               z             = z.new,
               n             = n,
               class         = classification,
               loglikelihood = loglike.new,
               ll.seq        = loglike[1:(j-1)],
               parameters.number=noparams,
               AIC           = AIC,
               BIC           = BIC,
               parameter.trace=list(alpha1.trace=alpha1.trace,
                                    alpha2.trace=alpha2.trace,
                                    alpha3.trace=alpha3.trace,
                                    beta.trace  =beta.trace),

               iterations    = (j-1),
               call          = tempcall)
  return(result)
}


MBGC.Estep.function <- function(data, alpha, beta){
  n <- dim(data)[1]
  s <- rep(0,n)
  logs <- rep(0,n)
  logy1s <- rep(0,n)
  logy2s <- rep(0,n)
  for (i in seq_len(n)) {
    expected.latent.res <- expected.latent(c(data[i,1],data[i,2]),
                                           alpha=c(alpha[1],alpha[2],alpha[3]),
                                           beta=beta)
    s[i] <- expected.latent.res$Expected.s
    logs[i] <- expected.latent.res$Expected.logs
    logy1s[i] <- expected.latent.res$Expected.logy1s
    logy2s[i] <- expected.latent.res$Expected.logy2s
    if ( is.nan(s[i]) || is.na(s[i]) ){stop("warings3")}
    #print(i)
  }
  return(list("s"=s, "logs"=logs, "logy1s"=logy1s, "logy2s"=logy2s))
}


MBGC.data.loglikelihood <- function(data, alpha1, alpha2, alpha3, beta, p.z){
  n <- dim(data)[1]
  G <- length(p.z)
  den <- matrix(0, nrow=n, ncol=G)
  for (i in seq_len(n)) {
    for (g in seq_len(G)) {
      den[i,g] <- dbivgamma(c(data[i,1],data[i,2]) ,
                            alpha=c(alpha1[g],alpha2[g],alpha3[g]),
                            beta=beta[g],
                            log=TRUE)$value
      if (den[i,g]==-Inf){stop("warings1")}
    }
  }
  newden <- sweep(den, 2, log(p.z), FUN="+")
  denom  <- matrixStats::rowLogSumExps(newden)
  return(sum(denom))
}


# try0 <- MGE(data=cutdata, maxit=100, tol=1e-5,
#             start.value = c(10,11,10,.01),
#             Aitken = FALSE, verbose=TRUE)
# try0$loglikelihood
# try0$BIC
#
# try2 <- MBGC(data=cutdata, G=2,
#              maxit=100, tol=1e-5,
#             start = "random.start",
#             Aitken = FALSE, verbose=TRUE)
# try2 <- MBGC(data=cutdata, G=2,
#              maxit=100, tol=1e-5,
#              start = "given",
#              start.value = matrix(c(.5,.5,.5,1e-4,.5,1.5, 1.5, 1.5, 1e-3,.5), ncol=2, byrow=FALSE),
#              Aitken = FALSE, verbose=TRUE)
# try2$loglikelihood
# try2$BIC
# try2$alpha1
# try2$alpha2
# try2$alpha3
# try2$beta
# try2$pro
#
#
# try3 <- MBGC(data=cutdata, G=3, maxit=100, tol=1e-5,
#              start = "random.start",
#              Aitken = FALSE, verbose=TRUE)
# try3 <- MBGC(data=ndat, G=3, maxit=100, tol=1e-5,
#              start = "given",
#              start.value = matrix(c(2,.5,2,1e-4,1/3,.5,2,2,1e-3,1/3,.5,.5,.5,1e-3,1/3), ncol=3, byrow=FALSE),
#              Aitken = FALSE, verbose=TRUE)
# try3$BIC
# try3$loglikelihood
# table(try3$class)
# plot(try3$loglike.seq, pch=20)
# try3$alpha1
# try3$alpha2
# try3$alpha3
# try3$beta
# try3$pro
#
# try4 <- MBGC(data=cutdata, G=4, maxit=100, tol=1e-5,
#              start = "random.start",
#              Aitken = FALSE, verbose=TRUE)
#
#
# plot(cutdata, pch=20, type="n")
# points(cutdata[which(try3$class==1),], pch=20)# , cex=(1-pzy[which(classification==1),1])*5)
# points(cutdata[which(try3$class==2),], pch=20, col=2)#, cex=(1-pzy[which(classification==2),2])*5)
# points(cutdata[which(try3$class==3),], pch=20, col=3)#, cex=(1-pzy[which(classification==2),2])*5)
#
#
# plot(cutdata, pch=20, type="n", xlab="AD", ylab="PD")
# points(cutdata[which(try2$class==1),], pch=20, col="brown4", cex=.6)
# points(cutdata[which(try2$class==2),], pch=20, col="darkgreen", cex=.6)
#
# ha <- Mclust(data=cutdata, G=2)
# plot(ha, what="classification")
# ha$loglik
# ha$bic
#
