#-------------------------------------
# Method 5: CVI / VVI / EVI 
#--------------------------------------

MBGR5 <- function(data,
                  G,
                  l1,
                  l2,
                  l3,
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
  } else {stop("l1 cannot not be NULL and has to be a formula")}
  if (!is.null(l2) && class(l2)=="formula"){
    if (as.character(l2[3])=="."){l2.new <- formula( paste( "x2", paste( names(data1),'',collapse='+',sep='' ), sep='~'))}
    else {l2.new <- formula(paste("x2", as.character(l2)[3], sep="~"))}
  } else {stop("l2 cannot not be NULL and has to be a formula")}
  if (!is.null(l3) && class(l3)=="formula"){
    if (as.character(l3[2])=="."){l3.new <- formula( paste( "x3", paste( names(data1), "", collapse = "+", sep = ""), sep="~"))}
    else {l3.new <- formula(paste("x3", as.character(l3)[2], sep="~"))}
  } else {stop("l3 cannot not be NULL and has to be a formula")}
  
  if (lp != "C" && lp != "E"){
    if (as.character(lp[2])=="."){lp.new <- formula( paste( "pzy", paste( names(data1), "", collapse = "+", sep = ""), sep="~"))}
    else {lp.new <- formula(paste("pzy", as.character(lp)[2], sep="~"))}
  }
  
  data1<-data; data1<-data1[ , names(data1)!=namey1]; data1<-data1[ , names(data1)!=namey2]
  if (lp != "C" && lp != "E") { datap<-data1 }
  
  #-----------------------------
  # starting values
  
  if (start=="mclust"){
    initial.mclust <- mclust::Mclust(G=G, data = cbind(y1,y2))
    z.current      <- initial.mclust$z
    
    if (lp == "C"){ 
      p.z.current        <- initial.mclust$parameters$pro 
      coef.p.current     <- NULL } else if (lp == "E"){
        p.z.current      <- rep(1/G, G) 
        coef.p.current   <- NULL } else {
          datap$pzy      <- z.current
          mp             <- nnet::multinom(formula = lp.new, data = datap, trace=FALSE)  # reltol
          coef.p.current <- coef(mp)
          p.z.current    <- stats::predict(mp, type="probs")   }
    
    index          <- initial.mclust$classification
    coef1.current  <- NULL
    coef2.current  <- NULL
    coef3.current  <- NULL
    alpha1.current <- matrix(0,ncol=G, nrow=n)
    alpha2.current <- matrix(0,ncol=G, nrow=n)
    alpha3.current <- matrix(0,ncol=G, nrow=n)
    
    Model.Matrix.1 <- as.matrix(Matrix::sparse.model.matrix(l1.new[-2], data=data1))
    Model.Matrix.2 <- as.matrix(Matrix::sparse.model.matrix(l2.new[-2], data=data1))
    Model.Matrix.3 <- as.matrix(Matrix::sparse.model.matrix(l3.new[-2], data=data1))
    
    x3 <- rep(0, n)
    for (i in seq_len(n)){ x3[i] <- runif(1, min=0, max=min(y1[i], y2[i])) }
    x1 <- y1 - x3
    x2 <- y2 - x3
    for (gg in seq_len(G)){
      data.start          <- data1[index==gg,]
      data.start$x1       <- x1[index==gg]
      data.start$x2       <- x2[index==gg]
      data.start$x3       <- x3[index==gg]
      start.temp.glm.1    <- glm(x1 ~ Model.Matrix.1[index==gg,-1], data=data.start, family=Gamma(link="log"))
      start.temp.glm.2    <- glm(x2 ~ Model.Matrix.2[index==gg,-1], data=data.start, family=Gamma(link="log"))
      start.temp.glm.3    <- glm(x3 ~ Model.Matrix.3[index==gg,-1], data=data.start, family=Gamma(link="log"))
      
      if (any(is.na(start.temp.glm.1$coefficient)==TRUE)){
        start.temp.glm.1.coef <- replace(as.vector(start.temp.glm.1$coefficient), which(is.na(start.temp.glm.1$coefficient)), 0)
      } else {start.temp.glm.1.coef <- start.temp.glm.1$coefficient}
      if (any(is.na(start.temp.glm.2$coefficient)==TRUE)){
        start.temp.glm.2.coef <- replace(as.vector(start.temp.glm.2$coefficient), which(is.na(start.temp.glm.2$coefficient)), 0)
      } else {start.temp.glm.2.coef <- start.temp.glm.2$coefficient}
      if (any(is.na(start.temp.glm.3$coefficient)==TRUE)){
        start.temp.glm.3.coef <- replace(as.vector(start.temp.glm.3$coefficient), which(is.na(start.temp.glm.3$coefficient)), 0)
      } else {start.temp.glm.3.coef <- start.temp.glm.3$coefficient}
      coef1.current       <- cbind(coef1.current, as.vector(start.temp.glm.1.coef))
      coef2.current       <- cbind(coef2.current, as.vector(start.temp.glm.2.coef))
      coef3.current       <- cbind(coef3.current, as.vector(start.temp.glm.3.coef))
      alpha1.current[,gg] <- rep(1/summary(start.temp.glm.1)$dispersion, n)
      alpha2.current[,gg] <- rep(1/summary(start.temp.glm.2)$dispersion, n)
      alpha3.current[,gg] <- rep(1/summary(start.temp.glm.3)$dispersion, n)
    }
    beta.current          <- sum( z.current*(alpha1.current+alpha2.current+alpha3.current) ) / sum( z.current * rep.col(x1+x2+x3, G) ) 
  }
  if (start=="random"){
    stop("start should be as Mclust (for now)")
  }
  
  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)
  loglike.store   <- rep(-.Machine$double.xmax, 3)
  
  coef1.trace     <- list(coef1.current)
  coef2.trace     <- list(coef2.current)
  coef3.trace     <- list(coef3.current)
  z.trace         <- list(z.current)
  coef.p.trace    <- coef.p.current
  
  alpha1.trace    <- list(alpha1.current)
  alpha2.trace    <- list(alpha2.current)
  alpha3.trace    <- list(alpha3.current)
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
                                                alpha=c(alpha1.current[i,g],alpha2.current[i,g],alpha3.current[i,g]),
                                                beta=beta.current)
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
    if (lp=="C" || lp=="E"){ 
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
      ait <- MoEClust::MoE_aitken(loglike.store)
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
      p.z.new        <- colSums(z.new)/n 
      coef.p.new     <- NULL  
      mp             <- NULL } else if (lp == "E"){
        p.z.new      <- rep(1/G, G) 
        coef.p.new   <- NULL  
        mp           <- NULL } else {
          datap$pzy  <- z.new
          mp         <- nnet::multinom(formula = lp.new, data = datap, trace=FALSE)  # reltol
          coef.p.new <- coef(mp)
          p.z.new    <- predict(mp, type="probs") }
    
    coef1.new  <- coef1.current
    coef2.new  <- coef2.current
    coef3.new  <- coef3.current
    alpha1.new <- matrix(0, ncol=G, nrow=n)
    alpha2.new <- matrix(0, ncol=G, nrow=n)
    alpha3.new <- matrix(0, ncol=G, nrow=n)
    
    hessian1   <- NULL
    hessian2   <- NULL
    hessian3   <- NULL
    
    Q.function <- function(coef, Model.Matrix, Expected.log.value, z.mat, wh){
      q.res <- sum( z.mat[,wh] * exp(coef%*%t(Model.Matrix)) * log(beta.current) ) - 
                sum( z.mat[,wh] * lgamma(exp(coef%*%t(Model.Matrix))) ) +
                sum( z.mat[,wh] * exp(coef%*%t(Model.Matrix)) * Expected.log.value[,wh] )
      return(q.res)
    }
    
    for (g in c(1:G)){
      m1             <- glm(alpha1.current[,g] ~ Model.Matrix.1[,-1], family=Gamma(link="log")  )#, weights=z.new[,g])
      coef1.new.temp <- as.vector(m1$coefficients)
      coef1.optim    <- optim(par               = coef1.new.temp,
                              fn                = Q.function, 
                              Model.Matrix      = Model.Matrix.1, 
                              Expected.log.value= Expected.logx1,
                              z.mat             = z.new,
                              wh                = g,
                              gr=NULL,  hessian = TRUE, 
                              control           = list(fnscale=-1, maxit=2e+9))
      if (coef1.optim$convergence != 0) cat("optim coef1 not converged at G=",g, ". \n")
      coef1.new[,g]  <- coef1.optim$par
      alpha1.new[,g] <- as.vector(exp(coef1.new[,g]%*%t(Model.Matrix.1)))  
      hessian1       <- append(hessian1, list(coef1.optim$hessian))
      
      m2             <- glm(alpha2.current[,g] ~ Model.Matrix.2[,-1], family=Gamma(link="log") )#, weights=z.new[,g]
      coef2.new.temp <- as.vector(m2$coefficients)
      coef2.optim    <- optim(par               = coef2.new.temp,
                              fn                = Q.function, 
                              Model.Matrix      = Model.Matrix.2, 
                              Expected.log.value= Expected.logx2,
                              z.mat             = z.new,
                              wh                = g,
                              gr=NULL,  hessian = TRUE, 
                              control           = list(fnscale=-1, maxit=2e+9))
      if (coef2.optim$convergence != 0) cat("optim coef2 not converged at G=",g, ". \n")
      coef2.new[,g]  <- coef2.optim$par
      alpha2.new[,g] <- as.vector(exp(coef2.new[,g]%*%t(Model.Matrix.2)))  
      hessian2       <- append(hessian2, list(coef2.optim$hessian))
      
      m3             <- glm(alpha3.current[,g] ~ Model.Matrix.3[,-1], family=Gamma(link="log") )#, weights=z.new[,g])
      coef3.new.temp <- as.vector(m3$coefficients)
      coef3.optim    <- optim(par               = coef3.new.temp,  
                              fn                = Q.function, 
                              Model.Matrix      = Model.Matrix.3, 
                              Expected.log.value= Expected.logx3,
                              z.mat             = z.new,
                              wh                = g,
                              gr=NULL,  hessian = TRUE, 
                              control           = list(fnscale=-1, maxit=2e+9))
      if (coef3.optim$convergence != 0) cat("optim coef3 not converged at G=",g, ". \n")
      coef3.new[,g]  <- coef3.optim$par
      alpha3.new[,g] <- as.vector(exp(coef3.new[,g]%*%t(Model.Matrix.3)))
      hessian3       <- append(hessian3, list(coef3.optim$hessian))
    }
    
    beta.new         <- sum( z.new*(alpha1.new + alpha2.new + alpha3.new) ) / sum( z.new*(Expected.x1+Expected.x2+Expected.x3) )
    
    if (verbose){
      cat(c("alpha1:", round(colMeans(alpha1.new),4), "\n"))
      cat(c("alpha2:", round(colMeans(alpha2.new),4), "\n"))
      cat(c("alpha3:", round(colMeans(alpha3.new),4), "\n"))
      cat(c("beta:",   round(beta.new,4), "\n"))
      if (is.matrix(p.z.new)) { cat(c("p.z:", round(colMeans(p.z.new),4), "\n")) } else cat(c("p.z:", round(p.z.new, 4), "\n"))
    }
    
    loglike.current<- loglike.new
    alpha1.current <- alpha1.new 
    alpha2.current <- alpha2.new 
    alpha3.current <- alpha3.new 
    beta.current   <- beta.new 
    p.z.current    <- p.z.new
    z.current      <- z.new
    coef1.current  <- coef1.new
    coef2.current  <- coef2.new
    coef3.current  <- coef3.new
    coef.p.current <- coef.p.new
    
    alpha1.trace   <- append(alpha1.trace, list(alpha1.new))
    alpha2.trace   <- append(alpha2.trace, list(alpha2.new))
    alpha3.trace   <- append(alpha3.trace, list(alpha3.new))
    beta.trace     <- c(beta.trace, beta.new)
    coef1.trace    <- append(coef1.trace, list(coef1.new))
    coef2.trace    <- append(coef2.trace, list(coef2.new))
    coef3.trace    <- append(coef3.trace, list(coef3.new))
    coef.p.trace   <- append(coef.p.trace, list(coef.p.new))
    
    j              <- j+1
  }
  
  if (!all(loglike[1:(j-1)] == cummax(loglike[1:(j-1)]))){
    cat("Loglikelihood is not strictly monotonically increasing!", "\n")
  }
  
  all.hessian <- c(hessian1, hessian2, hessian3)
  if (!all(sapply(all.hessian, matrixcalc::is.negative.definite))){
    cat("Hessian matrix is not negative-definite.", "\n")
  }
  
  #	calculation of BIC and AIC for bivpoisson model
  if (lp=="C"){
    noparams <- G*(dim(coef1.current)[1]) + G*(dim(coef2.current)[1]) + G*(dim(coef3.current)[1]) + 1 + (G-1)
  } else if (lp=="E"){
    noparams <- G*(dim(coef1.current)[1]) + G*(dim(coef2.current)[1]) + G*(dim(coef3.current)[1]) + 1
  } else {
    noparams <- G*(dim(coef1.current)[1]) + G*(dim(coef2.current)[1]) + G*(dim(coef3.current)[1]) + 1 + mp$rank 
  }
  
  AIC<- -2*loglike[j-1] + noparams * 2
  BIC<- -2*loglike[j-1] + noparams * log(n)
  
  classification <- apply(z.current, 1, which.max)
  
  # fitted values
  if (lp=="C" || lp=="E"){
    y1.fitted <- sweep((alpha1.current+alpha3.current),2,rep(beta.current,G),FUN="/") %*% p.z.current
    y2.fitted <- sweep((alpha2.current+alpha3.current),2,rep(beta.current,G),FUN="/") %*% p.z.current
  }  else {
    y1.fitted <- rowSums(p.z.current * sweep((alpha1.current+alpha3.current),2,rep(beta.current,G),FUN="/") )
    y2.fitted <- rowSums(p.z.current * sweep((alpha2.current+alpha3.current),2,rep(beta.current,G),FUN="/") )
  }
  y1.residual <- y1 - y1.fitted
  y2.residual <- y2 - y2.fitted
  
  # formula
  l1.formula <- formula(paste(namey1, as.character(l1.new)[3], sep="~"))
  l2.formula <- formula(paste(namey2, as.character(l2.new)[3], sep="~"))
  l3.formula <- formula(paste("", as.character(l2.new)[3], sep="~"))
  if (lp == "C") {
    lp.formula <- NULL } else if (lp == "E") {
      lp.formula <- "E"} else {
        lp.formula <- formula(paste("", as.character(lp.new)[3], sep="~"))}
  
  #  Calculation of output
  result<-list(coefficients = c(list(coef1.current), 
                                list(coef2.current), 
                                list(coef3.current),
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
               
               loglikelihood= loglike.current,
               llseq        = loglike[1:(j-1)],
               parameters.number=noparams,
               AIC          = AIC,
               BIC          = BIC,
               
               Hessian      = list(hessian1 = hessian1,
                                   hessian2 = hessian2,
                                   hessian3 = hessian3),
               
               parameter.trace=list(alpha1.trace=alpha1.trace,
                                    alpha2.trace=alpha2.trace,
                                    alpha3.trace=alpha3.trace,
                                    beta.trace  =beta.trace,  
                                    coef1.trace =coef1.trace, 
                                    coef2.trace =coef2.trace, 
                                    coef3.trace =coef3.trace, 
                                    coef.p.trace=coef.p.trace),
               
               n            = n,
               y            = cbind(y1, y2),
               Model.Matrix = list(Model.Matrix.1, Model.Matrix.2, Model.Matrix.3),
               formula      = list(l1.formula, l2.formula, l3.formula, lp.formula),
               p.model      = mp,
               
               call         = tempcall,
               iterations   = (j-1))
  options(warn=0)
  class(result)<-c('MBGR', 'glm')
  return(result)
}

