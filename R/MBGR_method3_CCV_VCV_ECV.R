#-------------------------------------
# Method 3: assume beta is regressed on its covariates, same alpha for all obs in each cluster
#--------------------------------------

MBGR3 <- function(data,
                  G,
                  l1,
                  l2,
                  l3,
                  l4,
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
                   l1     = l1,
                   l2     = l2,
                   l3     = l3,
                   l4     = l4,
                   lp     = lp,
                   expo   = expo,
                   maxit  = maxit,
                   tol    = tol,
                   Aitken = Aitken,
                   verbose= verbose)
  tempcall <- as.call( c(expression(MBGR), templist) ); rm(templist)
  #--------------------------------------------------
  
  namey1 <- as.character(l1[2])
  namey2 <- as.character(l2[2])
  y1     <- data[,names(data)==namey1]
  y2     <- data[,names(data)==namey2]
  n      <- length(y1)
  if (!is.null(expo)) {weight.var<-data[,names(data)==expo]} else {weight.var<-rep(1, n)}
  data$weight.var <- weight.var # there should be better ways of doing this
  
  if (!is.null(l1) && class(l1)=="formula"){
    if (as.character(l1[3])=="."){
      l1.new <- formula( paste( "x1", paste( names(data1),'',collapse='+',sep='' ), sep='~'))
    } else {l1.new <- formula(paste("x1", as.character(l1)[3], sep="~"))}
  }
  if (!is.null(l2) && class(l2)=="formula"){
    if (as.character(l2[3])=="."){l2.new <- formula( paste( "x2", paste( names(data1),'',collapse='+',sep='' ), sep='~'))}
    else {l2.new <- formula(paste("x2", as.character(l2)[3], sep="~"))}
  }
  if (!is.null(l3) && class(l3)=="formula"){
    if (as.character(l3[2])=="."){l3.new <- formula( paste( "x3", paste( names(data1), "", collapse = "+", sep = ""), sep="~"))}
    else {l3.new <- formula(paste("x3", as.character(l3)[2], sep="~"))}
  }
  if (!is.null(l4) && class(l4)=="formula"){
    if (as.character(l4[2])=="."){l4.new <- formula( paste( "be", paste( names(data1), "", collapse = "+", sep = ""), sep="~"))}
    else {l4.new <- formula(paste("be", as.character(l4)[2], sep="~"))}
  } else stop("l4 cannot be NULL and must be a formula")
  
  if (lp != "C" && lp != "E"){
    if (as.character(lp[2])=="."){lp.new <- formula( paste( "pzy", paste( names(data1), "", collapse = "+", sep = ""), sep="~"))}
    else {lp.new <- formula(paste("pzy", as.character(lp)[2], sep="~"))}
  }
  
  data1<-data; data1<-data1[ , names(data1)!=namey1]; data1<-data1[ , names(data1)!=namey2]
  if (lp != "C" && lp != "E") {datap<-data1}
  
  #-----------------------------
  # starting values
  
  if (start=="mclust"){
    initial.mclust <- mclust::Mclust(G=G, data = cbind(y1,y2))
    z.current      <- initial.mclust$z
    
    if (lp == "C"){ 
      p.z.current        <- initial.mclust$parameters$pro 
      coef.p.current     <- NULL } else if (lp == "E"){
        p.z.current      <- rep(1/G,G)
        coef.p.current   <- NULL } else {
          datap$pzy      <- z.current
          mp             <- nnet::multinom(formula = lp.new, data = datap, trace=FALSE)  # reltol
          coef.p.current <- coef(mp)
          p.z.current    <- stats::predict(mp, type="probs")   }

    index          <- initial.mclust$classification
    coef4.current  <- NULL
    alpha1.current <- rep(0,G)
    alpha2.current <- rep(0,G)
    alpha3.current <- rep(0,G)
    beta.current   <- matrix(0,ncol=G, nrow=n)
    
    Model.Matrix.4 <- as.matrix(sparse.model.matrix(l4.new[-2], data=data1))
    
    x3 <- rep(0, n)
    for (i in seq_len(n)){ x3[i] <- runif(1, min=0, max=min(y1[i], y2[i])) }
    x1 <- y1 - x3
    x2 <- y2 - x3
    for (gg in seq_len(G)){
      data.start          <- data1[index==gg,]
      data.start$x1       <- x1[index==gg]
      data.start$x2       <- x2[index==gg]
      data.start$x3       <- x3[index==gg]
      start.temp.glm.1    <- glm(l1.new, data=data.start, family=Gamma(link="log"))
      start.temp.glm.2    <- glm(l2.new, data=data.start, family=Gamma(link="log"))
      start.temp.glm.3    <- glm(l3.new, data=data.start, family=Gamma(link="log"))
      alpha1.current[gg]  <- 1/summary(start.temp.glm.1)$dispersion
      alpha2.current[gg]  <- 1/summary(start.temp.glm.2)$dispersion
      alpha3.current[gg]  <- 1/summary(start.temp.glm.3)$dispersion
      be                  <- (alpha1.current[gg]+alpha2.current[gg]+alpha3.current[gg]) / (x1+x2+x3)
      data1$be            <- be
      start.temp.glm.4    <- glm(l4.new, data=data1, family=Gamma(link="log"))
      beta.current[,gg]   <- start.temp.glm.4$fitted.values
      coef4.current       <- cbind(coef4.current, as.vector(start.temp.glm.4$coefficients))
      #beta.current[,gg]   <- exp( as.vector(start.temp.glm.4$coefficients) %*% t(Model.Matrix.4) )
    }
  }
  if (start=="random") stop("start should be as Mclust")

  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)
  loglike.store   <- rep(-.Machine$double.xmax, 3)
  
  coef4.trace  <- list(coef4.current)
  coef.p.trace <- list(coef.p.current)
  z.trace      <- list(z.current)
  
  alpha1.trace <- alpha1.current
  alpha2.trace <- alpha2.current
  alpha3.trace <- alpha3.current
  beta.trace   <- list(beta.current)
  
  j            <- 1
  
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
                                                beta=beta.current[i,g])
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
    if (lp == "C" || lp == "E"){ 
      newden   <- sweep(den, 2, log(p.z.current), FUN="+") } else{
        newden <- den + log(p.z.current)  }
    denom   <- matrixStats::rowLogSumExps(newden)
    z.new   <- exp(sweep(newden, 1, denom, FUN="-"))
    z.trace <- append(z.trace, list(z.new))
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
    if (lp == "C"){ 
      p.z.new    <- colSums(z.new)/n 
      coef.p.new <- NULL
      mp         <- NULL } else if (lp=="E"){
        p.z.new    <- rep(1/G, G)
        coef.p.new <- NULL
        mp         <- NULL } else {
          datap$pzy  <- z.new
          mp         <- nnet::multinom(formula = lp.new, data = datap, trace=FALSE)  # reltol
          coef.p.new <- coef(mp)
          p.z.new    <- stats::predict(mp, type="probs") }
    
    coef4.new  <- coef4.current
    alpha1.new <- rep(0, G)
    alpha2.new <- rep(0, G)
    alpha3.new <- rep(0, G)
    beta.new   <- matrix(0, ncol=G, nrow=n)
    
    hessian4   <- NULL
    
    Q.function    <- function(coef, wh){
      q.b.res <- (alpha1.current[wh] + alpha2.current[wh] + alpha3.current[wh]) * sum( z.new[,wh] * (coef %*% t(Model.Matrix.4))) - 
        sum( z.new[,wh] * exp(coef %*% t(Model.Matrix.4)) * (Expected.x1[,wh] + Expected.x2[,wh] + Expected.x3[,wh]) )
      return(q.b.res)
    }
    
    for (g in c(1:G)){
      m4             <- glm(beta.current[,g] ~ Model.Matrix.4[,-1], family=Gamma(link="log")  )#, weights=z.new[,g])
      coef4.new.temp <- as.vector(m4$coefficients)
      coef4.optim    <- optim(par               = coef4.new.temp, #coef4.current[,g], 
                              fn                = Q.function,
                              wh                = g,
                              gr=NULL,  hessian = TRUE, 
                              control           = list(fnscale=-1, maxit=2e+9))
      if (coef4.optim$convergence != 0) cat("optim coef4 not converged at G=",g, ". \n")
      coef4.new[,g]  <- coef4.optim$par
      beta.new[,g]   <- as.vector(exp(coef4.new[,g]%*%t(Model.Matrix.4)))  
      hessian4       <- append(hessian4, list(coef4.optim$hessian))
      
      alpha1.rootfun <- function(alpha1.var, wh){
        z.new[,wh]%*%log(beta.new[,wh]) - digamma(alpha1.var)*sum(z.new[,wh]) + z.new[,wh]%*%Expected.logx1[,wh]
      }
      alpha1.new[g]  <- uniroot(f    = alpha1.rootfun,
                                wh   = g,
                                lower= sqrt(.Machine$double.eps),
                                upper= 100000,
                                tol  = sqrt(.Machine$double.xmin))$root
      alpha2.rootfun <- function(alpha2.var, wh){
        z.new[,wh]%*%log(beta.new[,wh]) - digamma(alpha2.var)*sum(z.new[,wh]) + z.new[,wh]%*%Expected.logx2[,wh]
      }
      alpha2.new[g]  <- uniroot(f    = alpha2.rootfun,
                                wh   = g,
                                lower= sqrt(.Machine$double.eps),
                                upper= 100000,
                                tol  = sqrt(.Machine$double.xmin))$root
      alpha3.rootfun <- function(alpha3.var, wh){
        z.new[,wh]%*%log(beta.new[,wh]) - digamma(alpha3.var)*sum(z.new[,wh]) + z.new[,wh]%*%Expected.logx3[,wh]
      }
      alpha3.new[g]  <- uniroot(f    = alpha3.rootfun,
                                wh   = g,
                                lower= sqrt(.Machine$double.eps),
                                upper= 100000,
                                tol  = sqrt(.Machine$double.xmin))$root
    }
    
    if (verbose){
      cat(c("alpha1:", round(alpha1.new,4), "\n"))
      cat(c("alpha2:", round(alpha2.new,4), "\n"))
      cat(c("alpha3:", round(alpha3.new,4), "\n"))
      cat(c("beta:", round(colMeans(beta.new),4), "\n"))
      if (is.matrix(p.z.new)) { cat(c("p.z:", round(colMeans(p.z.new),4), "\n")) } else cat(c("p.z:", round(p.z.new,4), "\n"))
    }
    
    loglike.current<- loglike.new
    alpha1.current <- alpha1.new 
    alpha2.current <- alpha2.new 
    alpha3.current <- alpha3.new 
    beta.current   <- beta.new 
    p.z.current    <- p.z.new
    z.current      <- z.new
    coef4.current  <- coef4.new

    alpha1.trace   <- rbind(alpha1.trace, alpha1.new)
    alpha2.trace   <- rbind(alpha2.trace, alpha2.new)
    alpha3.trace   <- rbind(alpha3.trace, alpha3.new)
    beta.trace     <- append(beta.trace, list(beta.new))
    coef4.trace    <- append(coef4.trace, list(coef4.new))
    coef.p.trace   <- append(coef.p.trace, list(coef.p.new))
    
    j              <- j+1
  }
  
  if (!all(loglike[1:(j-1)] == cummax(loglike[1:(j-1)]))){
    cat("Loglikelihood is not strictly monotonically increasing!", "\n")
  }
  
  all.hessian <- c(hessian4)
  if (!all(sapply(all.hessian, matrixcalc::is.negative.definite))){
    cat("Hessian matrix is not negative-definite.", "\n")
  }
  
  #	calculation of BIC and AIC for bivpoisson model
  if (lp == "C"){
    noparams<- G*(dim(coef4.current)[1]) + 3*G + (G-1)
  } else if (lp == "E"){
    noparams<- G*(dim(coef4.current)[1]) + 3*G
  } else {
    noparams<- G*(dim(coef4.current)[1]) + 3*G + mp$rank 
  }
  
  AIC <- -2*loglike[j-1] + noparams * 2
  BIC <- -2*loglike[j-1] + noparams * log(n)
  
  classification <- apply(z.current, 1, which.max)
  
  # fitted values
  if (lp == "C" || lp == "E"){
    y1.fitted <- sweep((1/beta.current), 2, (alpha1.current+alpha3.current), FUN="*") %*% p.z.current 
    y2.fitted <- sweep((1/beta.current), 2, (alpha2.current+alpha3.current), FUN="*") %*% p.z.current 
  }  else {
    y1.fitted <- rowSums(p.z.current * sweep((1/beta.current), 2, (alpha1.current+alpha3.current), FUN="*") )
    y2.fitted <- rowSums(p.z.current * sweep((1/beta.current), 2, (alpha2.current+alpha3.current), FUN="*") )
  }
  y1.residual <- y1 - y1.fitted
  y2.residual <- y2 - y2.fitted
  
  # formula
  l1.formula <- formula(paste(namey1, as.character(l1.new)[3], sep="~"))
  l2.formula <- formula(paste(namey2, as.character(l2.new)[3], sep="~"))
  l3.formula <- formula(paste("", as.character(l3.new)[3], sep="~"))
  l4.formula <- formula(paste("", as.character(l4.new)[3], sep="~"))
  if (lp == "C") {
    lp.formula <- "C" } else if (lp == "E") {
      lp.formula <- "E"} else {
        lp.formula <- formula(paste("", as.character(lp.new)[3], sep="~"))} 
  
  #  Calculation of output
  result<-list(coefficients = c(list(coef4.current), 
                                list(coef.p.current)),
               alpha1       = alpha1.current,
               alpha2       = alpha2.current,
               alpha3       = alpha3.current,
               beta         = beta.current,
               
               clusterProb  = p.z.current,
               z            = z.current,
               class        = classification,
               
               fitted.values= data.frame(y1=y1.fitted,
                                         y2=y2.fitted),
               
               loglikelihood= loglike[j-1],
               llseq        = loglike[1:(j-1)],
               parameters.number=noparams,
               AIC          = AIC,
               BIC          = BIC,
               
               Hessian      = list(hessian4 = hessian4),
               
               n            = n,
               y            = cbind(y1, y2),
               Model.Matrix = Model.Matrix.4 ,
               formula      = list(l1.formula, l2.formula, l3.formula, l4.formula, lp.formula),
               p.model      = mp,
               
               call         = tempcall,
               iterations   = (j-1))
  options(warn=0)
  class(result)<-c('MBGR', 'glm')
  return(result)
}

