expected.latent.predict <- function(y, alpha, beta){
  y1 <- y[1]
  y2 <- y[2]
  alpha1 <- alpha[1]
  alpha2 <- alpha[2]
  alpha3 <- alpha[3]
  beta <- beta
  
  #-------------------------
  int.numerator.s.res <- log(alpha3)-log(beta)+dbivgamma(y, c(alpha1, alpha2, alpha3+1), beta, log = TRUE)$value
  #---------------------------
  int.denominator.res <- dbivgamma(y, alpha, beta, log=TRUE)$value
  #------------------------------
  
  Expected.s      <- exp(int.numerator.s.res - int.denominator.res)
  
  if (is.nan(Expected.s)){stop("NaN value returned for Expected.s")}
  if (is.na(Expected.s)){stop("NA value returned for Expected.s")}
  
  return(list("Expected.s"=Expected.s, 
              "s.numerator"=int.numerator.s.res,
              "denominator"=int.denominator.res))
}



predict.BGR <- function(object,
                        newdata,
                        maxit=100,
                        tol=1e-5,
                        verbose=TRUE,
                        Aitken=FALSE){
  
  #options(warn=-1)
  #tempcall
  n      <- dim(newdata)[1]
  
  coef1 <- object$coefficients[[1]]
  coef2 <- object$coefficients[[2]]
  coef3 <- object$coefficients[[3]]
  
  l1 <- object$formula[[1]]
  l2 <- object$formula[[2]]
  l3 <- object$formula[[3]]
  
  namey1 <- as.character(l1[2])
  namey2 <- as.character(l2[2])
  y1     <-newdata[,names(newdata)==namey1]
  y2     <-newdata[,names(newdata)==namey2]
  
  mod.mat.1 <- as.matrix(sparse.model.matrix(l1[-2], data=newdata))
  mod.mat.2 <- as.matrix(sparse.model.matrix(l2[-2], data=newdata))
  mod.mat.3 <- as.matrix(sparse.model.matrix(l3, data=newdata))
  
  # fitted alpha for newdata
  a1 <- as.vector(exp(coef1 %*% t(mod.mat.1)))
  a2 <- as.vector(exp(coef2 %*% t(mod.mat.2)))
  a3 <- as.vector(exp(coef3 %*% t(mod.mat.3)))
  
  # starting value
  b.current <- rep(mean(object$beta), n)
  
  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)
  loglike.store   <- rep(-.Machine$double.xmax, 3)
  
  j            <- 1
  
  while ( (loglike.diff > tol) && (j <= maxit) ) {
    
    #-------
    #E step
    #-------
    Expected.x3    <- rep(0,n)
    Expected.x1    <- rep(0,n)
    Expected.x2    <- rep(0,n)
    den            <- rep(0,n)
    for (i in seq_len(n)){
      expected.latent.res <- expected.latent.predict(c(y1[i],y2[i]),
                                                     alpha=c(a1[i],a2[i],a3[i]),
                                                     beta=beta.current[i])
      Expected.x3[i]      <- expected.latent.res$Expected.s
      Expected.x1[i]      <- y1[i]-expected.latent.res$Expected.s
      Expected.x2[i]      <- y2[i]-expected.latent.res$Expected.s
      if (Expected.x1[i] <=0 || Expected.x2[i] <=0){
        Expected.x3[i]  <- min(y1[i], y2[i])/2
        Expected.x1[i]  <- y1[i]-Expected.x3[i]
        Expected.x2[i]  <- y2[i]-Expected.x3[i]
        warning("negative expected x1 or x2 at ", i, " obs at ", j, " iteration", "\n")
      }
      den[i]              <- expected.latent.res$denominator
    }
    
    loglike.new    <- sum(den)
    loglike[j]     <- loglike.new
    loglike.store  <- c(loglike.store[-1], loglike.new)
    if (Aitken){
      loglike.diff <- abs( MoEClust::MoE_aitken(loglike.store)$linf - loglike.current )
    } else {
      loglike.diff <- abs(loglike.new-loglike.current)/(1+abs(loglike.new))
    }
    if (verbose){
      cat(c("iter:", j, "; current loglike:", round(loglike.new,6), "; loglike change:", round(loglike.diff,6), "\n"))
    }
    
    #--------
    # M step
    #--------
    beta.new <- (a1 + a2 + a3) / (Expected.x1 + Expected.x2 + Expected.x3)
    
    
    if (verbose){
      cat(c("beta:", summary(beta.new), "\n"))
    }
    
    loglike.current<- loglike.new
    beta.current   <- beta.new 
    #beta.trace     <- cbind(beta.trace, beta.new)
    j              <- j+1
  }
  
  newdata.y1.fitted <-  ((a1 + a3) / beta.current)
  newdata.y2.fitted <-  ((a2 + a3) / beta.current)
  
  result <-list(newdata.fit = data.frame(y1=newdata.y1.fitted, y2=newdata.y2.fitted),
                beta        = beta.current,
                alpha1      = a1,
                alpha2      = a2,
                alpha3      = a3,
                #beta.trace  =beta.trace,
                loglikelihood=loglike[j-1],
                llseq       = loglike[1:(j-1)],
                #call=tempcall,
                iterations=j-1)
  return(result)
}

pred.try <- predict.BGR(object=aviva1.train.copy, newdata=test.data)
plot(test.data$COST_AD, test.data$COST_TPD)
points(pred.try$newdata.fit, pch=20, col=2)

dim(test.data)

sim.y1 <- rep(0, 624)
sim.y2 <- rep(0, 624)
for (i in seq_len(624)){
  sim.y1[i] <- rgamma(1, shape=pred.try$alpha1[i]+pred.try$alpha3[i], pred.try$beta[i])
  sim.y2[i] <- rgamma(1, shape=pred.try$alpha2[i]+pred.try$alpha3[i], pred.try$beta[i])
}
plot(test.data$COST_AD, test.data$COST_TPD)
plot(test.data$COST_AD, test.data$COST_TPD)
points(sim.y1, sim.y2, pch=20, col=2)

#---------------------------------------------------------------------

aviva1.train
coef1 <- aviva1.train$coefficients[[1]]
coef2 <- aviva1.train$coefficients[[2]]
coef3 <- aviva1.train$coefficients[[3]]

MBGR.predict <- function(coef1, coef2, coef3, test.data){
  
  y1 <- test.data$COST_AD 
  y2 <- test.data$COST_TPD
  n <- dim(test.data)[1]
  
  mod.mat.1 <- as.matrix(sparse.model.matrix(~PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor + PH_PENPTS + MD_LIC_CAT_A, data=test.data))
  mod.mat.2 <- as.matrix(sparse.model.matrix(~PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor + PH_PENPTS + MD_LIC_CAT_A, data=test.data))
  mod.mat.3 <- as.matrix(sparse.model.matrix(~PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor + PH_PENPTS + MD_LIC_CAT_A, data=test.data))
  
  a1 <- as.vector(exp(coef1 %*% t(mod.mat.1)))
  a2 <- as.vector(exp(coef2 %*% t(mod.mat.2)))
  a3 <- as.vector(exp(coef3 %*% t(mod.mat.3)))
  
  # starting value
  b.current <- rep(mean(aviva1.train$beta), n)
  
  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)
  loglike.store   <- rep(-.Machine$double.xmax, 3)
  
  beta.trace   <- beta.current
  
  j            <- 1
  
  while ( (loglike.diff > tol) && (j <= maxit) ) {
    
    #-------
    #E step
    #-------
    Expected.x3    <- rep(0,n)
    Expected.x1    <- rep(0,n)
    Expected.x2    <- rep(0,n)
    den            <- rep(0,n)
    for (i in seq_len(n)){
      expected.latent.res <- expected.latent.predict(c(y1[i],y2[i]),
                                             alpha=c(a1[i],a2[i],a3[i]),
                                             beta=beta.current[i])
      Expected.x3[i]      <- expected.latent.res$Expected.s
      Expected.x1[i]      <- y1[i]-expected.latent.res$Expected.s
      Expected.x2[i]      <- y2[i]-expected.latent.res$Expected.s
      if (Expected.x1[i] <=0 || Expected.x2[i] <=0){
        Expected.x3[i]  <- min(y1[i], y2[i])/2
        Expected.x1[i]  <- y1[i]-Expected.x3[i]
        Expected.x2[i]  <- y2[i]-Expected.x3[i]
        warning("negative expected x1 or x2 at ", i, " obs at ", j, " iteration", "\n")
      }
      den[i]              <- expected.latent.res$denominator
    }
    
    loglike.new    <- sum(den)
    loglike[j]     <- loglike.new
    loglike.store  <- c(loglike.store[-1], loglike.new)
    if (Aitken){
      loglike.diff <- abs( MoEClust::MoE_aitken(loglike.store)$linf - loglike.current )
    } else {
      loglike.diff <- abs(loglike.new-loglike.current)/(1+abs(loglike.new))
    }
    if (verbose){
      cat(c("iter:", j, "; current loglike:", round(loglike.new,6), "; loglike change:", round(loglike.diff,6), "\n"))
    }
    
    #--------
    # M step
    #--------
    beta.new <- (a1 + a2 + a3) / (Expected.x1 + Expected.x2 + Expected.x3)
    
    
    if (verbose){
      cat(c("beta:", summary(beta.new), "\n"))
    }
    
    loglike.current<- loglike.new
    
    beta.current   <- beta.new 
    
    beta.trace     <- cbind(beta.trace, beta.new)
    
    j              <- j+1
    
  }
  
summary(beta.current)
y1.fitted <-  ((a1 + a3) / beta.current)
y2.fitted <-  ((a2 + a3) / beta.current)
}

plot(test.data$COST_AD, test.data$COST_TPD, main="mvGamma predictions", xlab="AD", ylab="PD")
points(y1.fitted, y2.fitted, pch=20, col=2)
plot(test.data$COST_AD, test.data$COST_TPD, main="univariate predictions", xlab="AD", ylab="PD")
points(pred1, pred2, pch=20, col=2)


mo1 <- glm(COST_AD ~ PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor + PH_PENPTS + MD_LIC_CAT_A, 
                    family = Gamma(link="log"), data=naviva.train)
mo2 <- glm(COST_TPD ~ PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor + PH_PENPTS + MD_LIC_CAT_A, 
                    family = Gamma(link="log"), data=naviva.train)
pred1 <- predict(mo1, newdata=test.data,  type="response", dispersion=1/MASS::gamma.shape(mo1))
pred2 <- predict(mo2, newdata=test.data,  type="response", dispersion=1/MASS::gamma.shape(mo2))

plot(cbind(naviva$COST_AD, naviva$COST_TPD), pch=20, cex=.7, xlab="AD", ylab="PD", col="cornflowerblue")
points(aviva.res$fitted.values, pch=20)
plot(cbind(naviva$COST_AD, naviva$COST_TPD), pch=20, cex=.7, xlab="AD", ylab="PD", col="cornflowerblue")
points(pred1, pred2, pch=20)

round(ComparisonMetrics(naviva.test$COST_AD, y1.fitted, expo = rep(1, 624)), 4)
round(ComparisonMetrics(naviva.test$COST_AD, pred1, expo = rep(1, 624)),4)

round(ComparisonMetrics(naviva.test$COST_TPD, y2.fitted, expo = rep(1, 624)), 4)
round(ComparisonMetrics(naviva.test$COST_TPD, pred2, expo = rep(1, 624)),4)

