#-------------------------------------
# Method 7: CVE / EVE / VVE
#--------------------------------------

MBGR.VE <- function(y,
                    data,
                    G,
                    f1,
                    f2,
                    f3,
                    f4,
                    gating,
                    initialization = "mclust",
                    maxit   = 200,
                    tol     = 1e-5,
                    verbose = TRUE){
  options(warn=-1)
  tempcall<-as.call(c(expression(MBGR),list(y      = y,
                                            data   = substitute(data),
                                            G      = G,
                                            f1     = f1,
                                            f2     = f2,
                                            f3     = f3,
                                            f4     = f4,
                                            gating = gating,
                                            initialization = initialization,
                                            maxit  = maxit,
                                            tol    = tol,
                                            verbose= verbose)))
  namey1 <- y[1]
  namey2 <- y[2]
  y1     <- data[,names(data)==namey1]
  y2     <- data[,names(data)==namey2]
  n      <- length(y1)

  newdata <- data
  newdata <- newdata[,names(newdata)!=namey1]
  newdata <- newdata[,names(newdata)!=namey2]

  if (!is.null(f1) && class(f1)=="formula"){
    if (as.character(f1[2])=="."){
      f1.new <- formula( paste( "x1", paste( names(newdata),'',collapse='+',sep='' ), sep='~'))
    } else {f1.new <- formula(paste("x1", as.character(f1)[2], sep="~"))}
  } else {stop("f1 cannot be NULL or has to be a formula.")}
  if (!is.null(f2) && class(f2)=="formula"){
    if (as.character(f2[2])=="."){f2.new <- formula( paste( "x2", paste( names(newdata),'',collapse='+',sep='' ), sep='~'))}
    else {f2.new <- formula(paste("x2", as.character(f2)[2], sep="~"))}
  } else {stop("f2 cannot be NULL or has to be a formula.")}
  if (!is.null(f3) && class(f3)=="formula"){
    if (as.character(f3[2])=="."){f3.new <- formula( paste( "x3", paste( names(newdata), "", collapse = "+", sep = ""), sep="~"))}
    else {f3.new <- formula(paste("x3", as.character(f3)[2], sep="~"))}
  } else {stop("f3 cannot be NULL or has to be a formula.")}
  if (!is.null(f4) && class(f4)=="formula"){
    if (as.character(f4[2])=="."){f4.new <- formula( paste( "be", paste( names(newdata), "", collapse = "+", sep = ""), sep="~"))}
    else {f4.new <- formula(paste("be", as.character(f4)[2], sep="~"))}
  } else {stop("f4 cannot be NULL or has to be a formula.")}
  if (gating!="C" && gating!="E"){
    if (as.character(gating[2])=="."){gating.new <- formula( paste( "pzy", paste( names(newdata), "", collapse = "+", sep = ""), sep="~"))}
    else {gating.new <- formula(paste("pzy", as.character(gating)[2], sep="~"))}
  }

  Model.Matrix.1 <- as.matrix(Matrix::sparse.model.matrix(f1.new[-2], data=newdata))
  Model.Matrix.2 <- as.matrix(Matrix::sparse.model.matrix(f2.new[-2], data=newdata))
  Model.Matrix.3 <- as.matrix(Matrix::sparse.model.matrix(f3.new[-2], data=newdata))
  Model.Matrix.4 <- as.matrix(Matrix::sparse.model.matrix(f4.new[-2], data=newdata))
  #-----------------------------
  # starting values

  if (initialization=="mclust"){
    initial.mclust <- mclust::Mclust(G=G, data = cbind(y1,y2), verbose = FALSE)
    z.init      <- initial.mclust$z

    if (gating=="C"){
      p.z.current        <- initial.mclust$parameters$pro
      coef.p.current     <- NULL } else if (gating=="E"){
        p.z.current      <- rep(1/G, G)
        coef.p.current   <- NULL } else {
          newdata$pzy      <- z.init
          mp             <- nnet::multinom(formula = gating.new, data = newdata, trace=FALSE)  # reltol
          coef.p.current <- coef(mp)
          p.z.current    <- stats::predict(mp, type="probs")   }

    index          <- initial.mclust$classification
    coef1.current  <- NULL
    coef2.current  <- NULL
    coef3.current  <- NULL
    coef4.current  <- NULL
    alpha1.current <- matrix(0,ncol=G, nrow=n)
    alpha2.current <- matrix(0,ncol=G, nrow=n)
    alpha3.current <- matrix(0,ncol=G, nrow=n)
    beta.current   <- rep(0,n)

    x3 <- rep(0, n)
    for (i in seq_len(n)){ x3[i] <- runif(1, min=0, max=min(y1[i], y2[i])) }
    x1 <- y1 - x3
    x2 <- y2 - x3
    for (gg in seq_len(G)){
      data.start          <- newdata[index==gg,]
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
    newdata$be          <- rowSums(( z.init*(alpha1.current+alpha2.current+alpha3.current) ) / ( z.init*rep.col(x1+x2+x3, G) ))
    start.temp.glm.4    <- glm(f4.new, data=newdata, family=Gamma(link="log"))
    beta.current        <- start.temp.glm.4$fitted.values
    coef4.current       <- as.vector(start.temp.glm.4$coefficients)
  }

  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)
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
                                                beta=beta.current[i])
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
      }
    }
    if (gating=="C" || gating=="E"){
      newden   <- sweep(den, 2, log(p.z.current), FUN="+") } else{
        newden <- den + log(p.z.current)  }
    denom   <- matrixStats::rowLogSumExps(newden)
    z.new   <- exp(sweep(newden, 1, denom, FUN="-"))

    loglike.new    <- sum(denom)
    loglike[j]     <- loglike.new
    loglike.diff <- abs(loglike.new-loglike.current)/(1+abs(loglike.new))
    if (verbose){
      cat(c("iter:", j, "; current loglike:", round(loglike.new,6), "; loglike change:", round(loglike.diff,6), "\n"))
    }

    #--------
    # M step
    #--------
    if (gating == "C"){
      p.z.new        <- colSums(z.new)/n
      coef.p.new     <- NULL
      mp             <- NULL } else if (gating == "E"){
        p.z.new      <- rep(1/G, G)
        coef.p.new   <- NULL
        mp           <- NULL } else {
          newdata$pzy  <- z.new
          mp         <- nnet::multinom(formula = gating.new, data = newdata, trace=FALSE)  # reltol
          coef.p.new <- coef(mp)
          p.z.new    <- stats::predict(mp, type="probs") }

    coef1.new  <- coef1.current
    coef2.new  <- coef2.current
    coef3.new  <- coef3.current
    alpha1.new <- matrix(0, ncol=G, nrow=n)
    alpha2.new <- matrix(0, ncol=G, nrow=n)
    alpha3.new <- matrix(0, ncol=G, nrow=n)
    beta.new   <- rep(0, n)
    hessian1   <- NULL
    hessian2   <- NULL
    hessian3   <- NULL
    hessian4   <- NULL
    Q.a.function <- function(coef, Model.Matrix, Expected.log.value, wh){
      q.a.res <- sum( z.new[,wh] * exp(coef%*%t(Model.Matrix)) * log(beta.current) ) -
                  sum( z.new[,wh] * lgamma(exp(coef%*%t(Model.Matrix))) ) +
                  sum( z.new[,wh] * exp(coef%*%t(Model.Matrix)) * Expected.log.value[,wh] )
      return(q.a.res)
    }
    Q.b.function <- function(coef){
      q.b.res <- sum( rowSums(z.new * (alpha1.new + alpha2.new + alpha3.new)) * (coef %*% t(Model.Matrix.4)) ) -
                   sum( rowSums(z.new * (Expected.x1 + Expected.x2 + Expected.x3)) * exp(coef %*% t(Model.Matrix.4)) )
      return(q.b.res)
    }
    for (g in c(1:G)){
      m1             <- glm(alpha1.current[,g] ~ Model.Matrix.1[,-1], family=Gamma(link="log"))
      coef1.optim    <- optim(par               = as.vector(m1$coefficients),
                              fn                = Q.a.function,
                              Model.Matrix      = Model.Matrix.1,
                              Expected.log.value= Expected.logx1,
                              wh                = g,
                              gr                = NULL,
                              hessian           = TRUE,
                              control           = list(fnscale=-1, maxit=2e+9))
      if (coef1.optim$convergence != 0) cat("optim coef1 not converged at G=",g, ". \n")
      coef1.new[,g]  <- coef1.optim$par
      alpha1.new[,g] <- as.vector(exp(coef1.new[,g]%*%t(Model.Matrix.1)))
      hessian1       <- append(hessian1, list(coef1.optim$hessian))

      m2             <- glm(alpha2.current[,g] ~ Model.Matrix.2[,-1], family=Gamma(link="log"))
      coef2.new.temp <- as.vector(m2$coefficients)
      coef2.optim    <- optim(par               = coef2.new.temp,
                              fn                = Q.a.function,
                              Model.Matrix      = Model.Matrix.2,
                              Expected.log.value= Expected.logx2,
                              wh                = g,
                              gr                = NULL,
                              hessian           = TRUE,
                              control           = list(fnscale=-1, maxit=2e+9))
      if (coef2.optim$convergence != 0) cat("optim coef2 not converged at G=",g, ". \n")
      coef2.new[,g]  <- coef2.optim$par
      alpha2.new[,g] <- as.vector(exp(coef2.new[,g]%*%t(Model.Matrix.2)))
      hessian2       <- append(hessian2, list(coef2.optim$hessian))

      m3             <- glm(alpha3.current[,g] ~ Model.Matrix.3[,-1], family=Gamma(link="log") )#, weights=z.new[,g])
      coef3.new.temp <- as.vector(m3$coefficients)
      coef3.optim    <- optim(par               = coef3.new.temp,
                              fn                = Q.a.function,
                              Model.Matrix      = Model.Matrix.3,
                              Expected.log.value= Expected.logx3,
                              wh                = g,
                              gr                = NULL,
                              hessian           = TRUE,
                              control           = list(fnscale=-1, maxit=2e+9))
      if (coef3.optim$convergence != 0) cat("optim coef3 not converged at G=",g, ". \n")
      coef3.new[,g]  <- coef3.optim$par
      alpha3.new[,g] <- as.vector(exp(coef3.new[,g]%*%t(Model.Matrix.3)))
      hessian3       <- append(hessian3, list(coef3.optim$hessian))
    }
    m4             <- glm(beta.current ~ Model.Matrix.4[,-1], family=Gamma(link="log")  )#, weights=z.new[,g])
    coef4.new.temp <- as.vector(m4$coefficients)
    coef4.optim    <- optim(par               = coef4.new.temp,
                            fn                = Q.b.function,
                            gr                = NULL,
                            hessian           = TRUE,
                            control           = list(fnscale=-1, maxit=2e+9))
    if (coef4.optim$convergence != 0) cat("optim coef4 not converged at G=",g, ". \n")
    coef4.new      <- coef4.optim$par
    beta.new       <- as.vector(exp(coef4.new %*% t(Model.Matrix.4)))
    hessian4       <- coef4.optim$hessian

    if (verbose){
      cat(c("alpha1:", round(colMeans(alpha1.new),4), "\n"))
      cat(c("alpha2:", round(colMeans(alpha2.new),4), "\n"))
      cat(c("alpha3:", round(colMeans(alpha3.new),4), "\n"))
      cat(c("beta:", round(mean(beta.new),4), "\n"))
      if (is.matrix(p.z.new)) { cat(c("p.z:", round(colMeans(p.z.new),4), "\n")) } else cat(c("p.z:", round(p.z.new,4), "\n"))
    }

    loglike.current<- loglike.new
    alpha1.current <- alpha1.new
    alpha2.current <- alpha2.new
    alpha3.current <- alpha3.new
    beta.current   <- beta.new
    p.z.current    <- p.z.new
    coef1.current  <- coef1.new
    coef2.current  <- coef2.new
    coef3.current  <- coef3.new
    coef4.current  <- coef4.new
    coef.p.current <- coef.p.new
    j              <- j+1
  }

  if (!all(loglike[1:(j-1)] == cummax(loglike[1:(j-1)]))){
    cat("Loglikelihood is not strictly monotonically increasing!", "\n")
  }
  all.hessian <- c(hessian1, hessian2, hessian3, list(hessian4))
  if (!all(sapply(all.hessian, matrixcalc::is.negative.definite))){
    cat("Hessian matrix is not negative-definite.", "\n")
  }

  #	calculation of BIC and AIC for bivpoisson model
  if (gating=="C"){
    noparams <- G*(dim(coef1.current)[1]) + G*(dim(coef2.current)[1]) + G*(dim(coef3.current)[1]) + length(coef4.current) + (G-1)
  } else if (gating=="E"){
    noparams <- G*(dim(coef1.current)[1]) + G*(dim(coef2.current)[1]) + G*(dim(coef3.current)[1]) + length(coef4.current)
  } else {
    noparams <- G*(dim(coef1.current)[1]) + G*(dim(coef2.current)[1]) + G*(dim(coef3.current)[1]) + length(coef4.current) + mp$rank
  }
  AIC<- -2*loglike[j-1] + noparams * 2
  BIC<- -2*loglike[j-1] + noparams * log(n)
  classification <- apply(z.new, 1, which.max)

  # fitted values
  if (gating=="C" || gating=="E"){
    y1.fitted <- ((alpha1.current + alpha3.current) / rep.col(beta.current, G)) %*% p.z.current
    y2.fitted <- ((alpha2.current + alpha3.current) / rep.col(beta.current, G)) %*% p.z.current
  }  else {
    y1.fitted <- rowSums(p.z.current * ((alpha1.current + alpha3.current)/ rep.col(beta.current, G) ) )
    y2.fitted <- rowSums(p.z.current * ((alpha2.current + alpha3.current)/ rep.col(beta.current, G) ) )
  }
  y1.residual <- y1 - y1.fitted
  y2.residual <- y2 - y2.fitted

  # formula
  f1.formula <- formula(paste(namey1, as.character(f1.new)[3], sep="~"))
  f2.formula <- formula(paste(namey2, as.character(f2.new)[3], sep="~"))
  f3.formula <- formula(paste("", as.character(f3.new)[3], sep="~"))
  f4.formula <- formula(paste("", as.character(f4.new)[3], sep="~"))
  if (gating == "C") {
    newgating <- gating.formula <- "C" } else if (gating == "E") {
      newgating <- gating.formula <- "E"} else {
        gating.formula <- formula(paste("", as.character(gating.new)[3], sep="~"))
        newgating <- "V"}

  result<-list(modelName    = "VE",
               gating       = newgating,
               coefficients = list(expert = list(alpha1=coef1.current,
                                                 alpha2=coef2.current,
                                                 alpha3=coef3.current,
                                                 beta  =coef4.current),
                                   gating = coef.p.current),
               alpha1       = alpha1.current,
               alpha2       = alpha2.current,
               alpha3       = alpha3.current,
               beta         = beta.current,
               pro          = p.z.current,
               z            = z.new,
               class        = classification,
               G            = G,
               fitted.values= data.frame(y1=y1.fitted,
                                         y2=y2.fitted),
               residuals    = cbind(y1.residual, y2.residual),
               loglike      = loglike.current,
               ll           = loglike[1:(j-1)],
               df           = noparams,
               AIC          = AIC,
               BIC          = BIC,
               Hessian      = list(hessian1 = hessian1,
                                   hessian2 = hessian2,
                                   hessian3 = hessian3,
                                   hessian4 = hessian4),
               n            = n,
               y            = cbind(y1, y2),
               Model.Matrix = list(Model.Matrix.1,
                                   Model.Matrix.2,
                                   Model.Matrix.3,
                                   Model.Matrix.4),
               formula      = list(f1.formula,
                                   f2.formula,
                                   f3.formula,
                                   f4.formula,
                                   gating.formula),
               gating.model = mp,
               call         = tempcall,
               iterations   = (j-1))
  options(warn=0)
  structure(result, class = c('MBGR'))
}
