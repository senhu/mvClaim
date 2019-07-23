#-------------------------------------
# Method 6: CIV / VIV / EIV
#--------------------------------------

MBGR_IV <- function(y,
                    data,
                    G,
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

  newdata<-data
  newdata<-newdata[,names(newdata)!=namey1]
  newdata<-newdata[,names(newdata)!=namey2]

  if (!is.null(f4) && class(f4)=="formula"){
    if (as.character(f4[2])=="."){f4.new <- formula( paste( "be", paste( names(newdata), "", collapse = "+", sep = ""), sep="~"))}
    else {f4.new <- formula(paste("be", as.character(f4)[2], sep="~"))}
  } else stop("f4 cannot be NULL and must be a formula")

  if (gating != "C" && gating != "E"){
    if (as.character(gating[2])=="."){gating.new <- formula( paste( "pzy", paste( names(newdata), "", collapse = "+", sep = ""), sep="~"))}
    else {gating.new <- formula(paste("pzy", as.character(gating)[2], sep="~"))}
  }

  Model.Matrix.4 <- as.matrix(stats::model.matrix(f4.new[-2], data=newdata))
  #-----------------------------
  # starting values

  if (initialization=="mclust"){
    initial.mclust <- mclust::Mclust(G=G, data = cbind(y1,y2), verbose = FALSE)
    z.init         <- initial.mclust$z

    if (gating == "C"){
      p.z.current        <- initial.mclust$parameters$pro
      coef.p.current     <- NULL } else if (gating == "E"){
        p.z.current      <- rep(1/G,G)
        coef.p.current   <- NULL } else {
          newdata$pzy      <- z.init
          mp             <- nnet::multinom(formula = gating.new, data = newdata, trace=FALSE)  # reltol
          coef.p.current <- coef(mp)
          p.z.current    <- stats::predict(mp, type="probs")   }

    index          <- initial.mclust$classification
    x3 <- rep(0, n)
    for (i in seq_len(n)){ x3[i] <- runif(1, min=0, max=min(y1[i], y2[i])) }
    x1 <- y1 - x3
    x2 <- y2 - x3
    coef4.current  <- NULL
    beta.current   <- matrix(0,ncol=G, nrow=n)
    alpha1.current.temp <- rep(0,G)
    alpha2.current.temp <- rep(0,G)
    alpha3.current.temp <- rep(0,G)
    for (gg in seq_len(G)){
      alpha1.current.temp[gg] <- MASS::fitdistr(x1[index==gg], "gamma")$estimate[1]
      alpha2.current.temp[gg] <- MASS::fitdistr(x2[index==gg], "gamma")$estimate[1]
      alpha3.current.temp[gg] <- MASS::fitdistr(x3[index==gg], "gamma")$estimate[1]
      be.temp             <- (alpha1.current.temp[gg]+alpha2.current.temp[gg]+alpha3.current.temp[gg]) / (x1+x2+x3)
      start.temp.glm.4    <- stats::glm(be.temp~Model.Matrix.4[,-1], family=Gamma(link="log"))
      beta.current[,gg]   <- start.temp.glm.4$fitted.values
      coef4.current       <- cbind(coef4.current, as.vector(start.temp.glm.4$coefficients))
    }
    alpha1.current        <- mean(alpha1.current.temp)
    alpha2.current        <- mean(alpha2.current.temp)
    alpha3.current        <- mean(alpha3.current.temp)
  }

  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)
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
        expected.latent.res  <- expected_latent(c(y1[i],y2[i]),
                                                alpha=c(alpha1.current,
                                                        alpha2.current,
                                                        alpha3.current),
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
      }
    }
    if (gating == "C" || gating == "E"){
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
      mp             <- NULL } else if (gating=="E"){
        p.z.new      <- rep(1/G, G)
        coef.p.new   <- NULL
        mp           <- NULL } else {
          newdata$pzy  <- z.new
          mp         <- nnet::multinom(formula = gating.new, data = newdata, trace=FALSE)  # reltol
          coef.p.new <- coef(mp)
          p.z.new    <- stats::predict(mp, type="probs") }

    coef4.new  <- coef4.current
    beta.new   <- matrix(0, ncol=G, nrow=n)
    hessian4   <- NULL
    Q.function    <- function(coef, wh){
      q.b.res <- (alpha1.current + alpha2.current + alpha3.current) * sum( z.new[,wh] * (coef %*% t(Model.Matrix.4))) -
        sum( z.new[,wh] * exp(coef %*% t(Model.Matrix.4)) * (Expected.x1[,wh] + Expected.x2[,wh] + Expected.x3[,wh]) )
      return(q.b.res)
    }
    for (g in c(1:G)){
      m4             <- stats::glm(beta.current[,g] ~ Model.Matrix.4[,-1], family=Gamma(link="log"))
      coef4.optim    <- stats::optim(par     = as.vector(m4$coefficients),
                              fn      = Q.function,
                              wh      = g,
                              gr      = NULL,
                              hessian = TRUE,
                              control = list(fnscale=-1, maxit=2e+9))
      if (coef4.optim$convergence != 0) cat("optim coef4 not converged at G=",g, ". \n")
      coef4.new[,g]  <- coef4.optim$par
      beta.new[,g]   <- as.vector(exp(coef4.new[,g]%*%t(Model.Matrix.4)))
      hessian4       <- append(hessian4, list(coef4.optim$hessian))
    }
    alpha1.rootfun <- function(alpha1.var){
      sum(z.new*log(beta.new)) - digamma(alpha1.var)*sum(z.new) + sum(z.new*Expected.logx1)
    }
    alpha1.new  <- stats::uniroot(f     = alpha1.rootfun,
                           lower = .Machine$double.eps, # sqrt(.Machine$double.eps)
                           upper = 100000,
                           tol   = sqrt(.Machine$double.xmin))$root
    alpha2.rootfun <- function(alpha2.var){
      sum(z.new*log(beta.new)) - digamma(alpha2.var)*sum(z.new) + sum(z.new*Expected.logx2)
    }
    alpha2.new <- stats::uniroot(f    = alpha2.rootfun,
                          lower= .Machine$double.eps,
                          upper= 100000,
                          tol  = sqrt(.Machine$double.xmin))$root
    alpha3.rootfun <- function(alpha3.var){
      sum(z.new*log(beta.new)) - digamma(alpha3.var)*sum(z.new) + sum(z.new*Expected.logx3)
    }
    alpha3.new <- stats::uniroot(f    = alpha3.rootfun,
                          lower= .Machine$double.eps,
                          upper= 100000,
                          tol  = sqrt(.Machine$double.xmin))$root

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
    coef4.current  <- coef4.new
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
  if (gating == "C"){
    noparams<- G*(dim(coef4.current)[1]) + 3 + (G-1)
  } else if (gating == "E"){
    noparams<- G*(dim(coef4.current)[1]) + 3
  } else {
    noparams<- G*(dim(coef4.current)[1]) + 3 + mp$rank
  }
  AIC <- -2*loglike[j-1] + noparams * 2
  BIC <- -2*loglike[j-1] + noparams * log(n)
  classification <- apply(z.new, 1, which.max)

  # fitted values
  if (gating == "C" || gating == "E"){
    y1.fitted <- sweep((1/beta.current), 2, (rep(alpha1.current+alpha3.current,G)), FUN="*") %*% p.z.current
    y2.fitted <- sweep((1/beta.current), 2, (rep(alpha2.current+alpha3.current,G)), FUN="*") %*% p.z.current
  }  else {
    y1.fitted <- rowSums(p.z.current * sweep((1/beta.current), 2, (rep(alpha1.current+alpha3.current,G)), FUN="*") )
    y2.fitted <- rowSums(p.z.current * sweep((1/beta.current), 2, (rep(alpha2.current+alpha3.current,G)), FUN="*") )
  }
  y1.residual <- y1 - y1.fitted
  y2.residual <- y2 - y2.fitted

  # formula
  f4.formula <- formula(paste("", as.character(f4.new)[3], sep="~"))
  if (gating == "C") {
    newgating <- gating.formula <- "C" } else if (gating == "E") {
      newgating <- gating.formula <- "E"} else {
        gating.formula <- formula(paste("", as.character(gating.new)[3], sep="~"))
        newgating <- "V"}

  #  Calculation of output
  result<-list(modelName    = "IV",
               gating       = newgating,
               coefficients = list(expert = list(coef4.current),
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
               loglike      = loglike[j-1],
               ll           = loglike[1:(j-1)],
               df           = noparams,
               AIC          = AIC,
               BIC          = BIC,
               Hessian      = list(hessian4),
               n            = n,
               y            = cbind(y1, y2),
               Model.Matrix = Model.Matrix.4 ,
               formula      = list(f4.formula, gating.formula),
               gating.model = mp,
               call         = tempcall,
               iterations   = (j-1))
  options(warn=0)
  structure(result, class = c('MBGR'))
}

