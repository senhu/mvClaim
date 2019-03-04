#-------------------------------------
# Gating only : VCC
#--------------------------------------

MBGR_Gating <- function(data,
                        response,
                        G,
                        lp,
                        expo    = NULL,
                        maxit   = 300,
                        tol     = 1e-6,
                        start   = "mclust",
                        Aitken  = FALSE,
                        verbose = TRUE){
  #-------------------------------------------------
  #options(warn=-1)
  templist <- list(data   = substitute(data),
                   G      = G,
                   response= response,
                   lp     = lp,
                   expo   = expo,
                   maxit  = maxit,
                   tol    = tol,
                   Aitken = Aitken,
                   verbose= verbose)
  tempcall <- as.call( c(expression(MBGR), templist) ); rm(templist)
  #--------------------------------------------------
  
  namey1 <- response[1]
  namey2 <- response[2]
  y1     <- data[,names(data)==namey1]
  y2     <- data[,names(data)==namey2]
  n      <- length(y1)
  
  if (!is.null(expo)) {weight.var<-data[,names(data)==expo]} else {weight.var<-rep(1, n)}
  data$weight.var <- weight.var # there should be better ways of doing this
  
  if (as.character(lp[2])=="."){lp.new <- formula( paste( "pzy", paste( names(data1), "", collapse = "+", sep = ""), sep="~"))}
    else {lp.new <- formula(paste("pzy", as.character(lp)[2], sep="~"))}

  data1 <- data; data1<-data1[ , names(data1)!=namey1]; data1<-data1[ , names(data1)!=namey2]
  datap <- data1
  
  #-----------------------------
  # starting values
  
  if (start=="mclust"){
    initial.mclust <- mclust::Mclust(G=G, data = cbind(y1,y2))
    z.current      <- initial.mclust$z
    
    datap$pzy      <- z.current
    mp             <- nnet::multinom(formula = lp.new, data = datap, trace=FALSE)  # reltol
    coef.p.current <- coef(mp)
    p.z.current    <- stats::predict(mp, type="probs")
    
    index          <- initial.mclust$classification
    alpha1.current <- rep(0,G)
    alpha2.current <- rep(0,G)
    alpha3.current <- rep(0,G)
    beta.current   <- rep(0,G)
    
    x3 <- rep(0, n)
    for (i in seq_len(n)){ x3[i] <- runif(1, min=0, max=min(y1[i], y2[i])) }
    x1 <- y1 - x3
    x2 <- y2 - x3
    for (gg in seq_len(G)){
      data.start          <- data1[index==gg,]
      data.start$x1       <- x1[index==gg]
      data.start$x2       <- x2[index==gg]
      data.start$x3       <- x3[index==gg]

      temp.v1 <- MASS::fitdistr(data.start$x1/mean(data.start$x1), "gamma")
      temp.v2 <- MASS::fitdistr(data.start$x2/mean(data.start$x2), "gamma")
      temp.v3 <- MASS::fitdistr(data.start$x3/mean(data.start$x3), "gamma")
      
      alpha1.current[gg] <- temp.v1$estimate[1]
      alpha2.current[gg] <- temp.v2$estimate[1]
      alpha3.current[gg] <- temp.v3$estimate[1]
      beta.current[gg] <- mean(c(temp.v1$estimate[2]/mean(data.start$x1),
                                 temp.v2$estimate[2]/mean(data.start$x2),
                                 temp.v3$estimate[2]/mean(data.start$x3)))
    }
  }
  if (start=="random"){
    stop("start should be as Mclust (for now)")
  }
  
  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)
  loglike.store   <- rep(-.Machine$double.xmax, 3)
  
  # z.trace         <- list(z.current)
  coef.p.trace    <- coef.p.current
  
  alpha1.trace    <- alpha1.current
  alpha2.trace    <- alpha2.current
  alpha3.trace    <- alpha3.current
  beta.trace      <- beta.current
  
  j               <- 1
  
  while ( (loglike.diff > tol) && (j <= maxit) ) {
    #-------
    #E step
    #-------
    Expected.x3    <- matrix(0,nrow=n, ncol=G)
    Expected.x1    <- matrix(0,nrow=n, ncol=G)
    Expected.x2    <- matrix(0,nrow=n, ncol=G)
    Expected.logx1 <- matrix(0,nrow=n, ncol=G)
    Expected.logx2 <- matrix(0,nrow=n, ncol=G)
    Expected.logx3 <- matrix(0,nrow=n, ncol=G)
    den            <- matrix(0,nrow=n, ncol=G)
    
    for (i in seq_len(n)){
      for (g in seq_len(G)){
        expected.latent.res  <- expected.latent(c(y1[i],y2[i]),
                                                alpha=c(alpha1.current[g],alpha2.current[g],alpha3.current[g]),
                                                beta=beta.current[g])
        Expected.x3[i,g]     <- expected.latent.res$Expected.s
        Expected.x1[i,g]     <- (y1[i]-expected.latent.res$Expected.s)
        Expected.x2[i,g]     <- (y2[i]-expected.latent.res$Expected.s)
        Expected.logx3[i,g]  <- expected.latent.res$Expected.logs
        Expected.logx1[i,g]  <- expected.latent.res$Expected.logy1s
        Expected.logx2[i,g]  <- expected.latent.res$Expected.logy2s
        if (Expected.x1[i,g] <=0 || Expected.x2[i,g] <=0){
          Expected.x3[i,g]  <- min(y1[i], y2[i])/2
          Expected.x1[i,g]  <- y1[i]-Expected.x3[i,g]
          Expected.x2[i,g]  <- y2[i]-Expected.x3[i,g]
          warning("negative expected x1 or x2 at ", i, " obs at ", j, " iteration", "\n")
        }
        den[i,g]             <- expected.latent.res$denominator
        #print(c(i,g))
      }
    }
    newden <- den + log(p.z.current)
    denom   <- matrixStats::rowLogSumExps(newden)
    z.new   <- exp(sweep(newden, 1, denom, FUN="-"))
    if (table(round(rowSums(z.new),0))!=n){stop("Wrong E-step estimation: z.new")}
    
    loglike.new    <- sum(denom)
    loglike[j]     <- loglike.new
    loglike.store  <- c(loglike.store[-1], loglike.new)
    if (Aitken){
      ait          <- MoEClust::MoE_aitken(loglike.store)
      loglike.diff <- ifelse(is.numeric(ait$a) && ait$a < 0, 0,abs(ait$linf-loglike.current))
      loglike.diff[is.nan(loglike.diff)] <- Inf
    } else {
      loglike.diff <- abs(loglike.new-loglike.current)/(1+abs(loglike.new))
    }
    if (verbose){
      cat(c("iter:", j, "; current loglike:", round(loglike.new,6), "; loglike change:", round(loglike.diff,6), "\n"))
    }
    
    
    #--------
    # M step
    #--------
    
    datap$pzy  <- z.new
    mp         <- nnet::multinom(formula = lp.new, data = datap, trace=FALSE)  # reltol
    coef.p.new <- coef(mp)
    p.z.new    <- stats::predict(mp, type="probs")
    
    alpha1.new <- rep(0, G)
    alpha2.new <- rep(0, G)
    alpha3.new <- rep(0, G)
    beta.new   <- rep(0, G)
    
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
    
    # hessian1   <- NULL
    # hessian2   <- NULL
    # hessian3   <- NULL
    
    if (verbose){
      cat(c("alpha1:", round(alpha1.new,4), "\n"))
      cat(c("alpha2:", round(alpha2.new,4), "\n"))
      cat(c("alpha3:", round(alpha3.new,4), "\n"))
      cat(c("beta:", round(beta.new,4), "\n"))
      cat(c("p.z:", round(colMeans(p.z.new),4), "\n"))
    }

    loglike.current<- loglike.new
    alpha1.current <- alpha1.new
    alpha2.current <- alpha2.new
    alpha3.current <- alpha3.new
    beta.current   <- beta.new
    p.z.current    <- p.z.new
    z.current      <- z.new
    coef.p.current <- coef.p.new
    
    alpha1.trace   <- rbind(alpha1.trace, alpha1.new)
    alpha2.trace   <- rbind(alpha2.trace, alpha2.new)
    alpha3.trace   <- rbind(alpha3.trace, alpha3.new)
    beta.trace     <- rbind(beta.trace, beta.new)
    coef.p.trace   <- append(coef.p.trace, list(coef.p.new))
    
    j              <- j+1
  }
  
  if (!all(loglike[1:(j-1)] == cummax(loglike[1:(j-1)]))){
    cat("Loglikelihood is not strictly monotonically increasing!", "\n")
  }
  
  # all.hessian <- c(hessian1, hessian2, hessian3)
  # if (!all(sapply(all.hessian, matrixcalc::is.negative.definite))){
  #   cat("Hessian matrix is not negative-definite.", "\n")
  # }

  #	calculation of BIC and AIC for bivpoisson model
  noparams <- G*4 + mp$rank

  AIC<- -2*loglike[j-1] + noparams * 2
  BIC<- -2*loglike[j-1] + noparams * log(n)
  
  classification <- apply(z.current, 1, which.max)
  
  # fitted values

  y1.fitted <- rowSums(sweep(p.z.current, 2, ((alpha1.current+alpha3.current)/beta.current), FUN="*"))
  y2.fitted <- rowSums(sweep(p.z.current, 2, ((alpha2.current+alpha3.current)/beta.current), FUN="*"))
  y1.residual <- y1 - y1.fitted
  y2.residual <- y2 - y2.fitted
  
  # formula
  lp.formula <- formula(paste("", as.character(lp.new)[3], sep="~"))
  
  #  Calculation of output
  result<-list(coefficients = list(coef.p.current),
               alpha1       = alpha1.current,
               alpha2       = alpha2.current,
               alpha3       = alpha3.current,
               beta         = beta.current,
               
               pro          = p.z.current,
               z            = z.current,
               class        = classification,
               
               fitted.values= data.frame(y1=y1.fitted,
                                         y2=y2.fitted),
               
               loglikelihood= loglike.current,
               llseq        = loglike[1:(j-1)],
               parameters.number=noparams,
               AIC          = AIC,
               BIC          = BIC,
               
               parameter.trace=list(alpha1.trace=alpha1.trace,
                                    alpha2.trace=alpha2.trace,
                                    alpha3.trace=alpha3.trace,
                                    beta.trace  =beta.trace,
                                    coef.p.trace=coef.p.trace),
               
               n            = n,
               y            = cbind(y1, y2),
               formula      = lp.formula,
               p.model      = mp,
               
               call         = tempcall,
               iterations   = (j-1))
  options(warn=0)
  class(result)<-c('MBGR', 'glm')
  return(result)
}

