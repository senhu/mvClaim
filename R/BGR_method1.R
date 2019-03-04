
BGR <- function(data,
                l1,
                l2,
                l3,
                expo=NULL,
                maxit=10,
                tol=1e-5,
                start,
                Aitken = TRUE,
                verbose=TRUE){
  
  #options(warn=-1)
  #tempcall
  namey1 <- as.character(l1[2])
  namey2 <- as.character(l2[2])
  y1     <- data[,names(data)==namey1]
  y2     <- data[,names(data)==namey2]
  n      <- length(y1)
  if (!is.null(expo)) {weight.var<-data[,names(data)==expo]} else {weight.var<-rep(1, n)}
  data$weight.var <- weight.var # there should be better ways of doing this
  
  if (!is.null(l1)){
    if (as.character(l1[3])=="."){
      l1.new <- formula( paste( "x1", paste( names(data1),'',collapse='+',sep='' ), sep='~')) 
    } else {l1.new <- formula(paste("x1", as.character(l1)[3], sep="~"))}
  } else stop("l1 cannot be NULL, must be supplied")
  if (!is.null(l2)){
    if (as.character(l2[3])=="."){l2.new <- formula( paste( "x2", paste( names(data1),'',collapse='+',sep='' ), sep='~'))}
    else {l2.new <- formula(paste("x2", as.character(l2)[3], sep="~"))}
  } else stop("l2 cannot be NULL, must be supplied")
  if (!is.null(l3)){
    if (as.character(l3[2])=="."){l3.new <- formula( paste( "x3", paste( names(data1), "", collapse = "+", sep = ""), sep="~"))}
    else {l3.new <- formula(paste("x3", as.character(l3)[2], sep="~"))}
  } else stop("l3 cannot be NULL, must be supplied")

  data1<-data; data1<-data1[ , names(data1)!=namey1]; data1<-data1[ , names(data1)!=namey2]
  datap<-data1
  
  #-----------------------------
  # starting values
  x3.start <- rep(0, n)
  #set.seed(10)
  for (i in seq_len(n)){
    x3.start[i] <- runif(1, min=0, max=min(y1[i], y2[i]))
  }
  x1.start       <- y1 - x3.start
  x2.start       <- y2 - x3.start
  data.start     <- data1
  data.start$x1  <- x1.start
  data.start$x2  <- x2.start
  data.start$x3  <- x3.start
  start.temp.glm.1 <- glm(l1.new, data=data.start, family=Gamma(link="log"))
  start.temp.glm.2 <- glm(l2.new, data=data.start, family=Gamma(link="log"))
  start.temp.glm.3 <- glm(l3.new, data=data.start, family=Gamma(link="log"))
  
  Model.Matrix.1 <- model.matrix(start.temp.glm.1)
  Model.Matrix.2 <- model.matrix(start.temp.glm.2)
  Model.Matrix.3 <- model.matrix(start.temp.glm.3)
  dispersion1.current <- summary(start.temp.glm.1)$dispersion
  dispersion2.current <- summary(start.temp.glm.2)$dispersion
  dispersion3.current <- summary(start.temp.glm.3)$dispersion
  coef1.current  <- start.temp.glm.1$coefficient
  coef2.current  <- start.temp.glm.2$coefficient
  coef3.current  <- start.temp.glm.3$coefficient
  
  alpha1.current <- rep(1/summary(start.temp.glm.1)$dispersion, n) #rep(MASS::gamma.shape(start.temp.glm.1)$alpha, n)
  alpha2.current <- rep(1/summary(start.temp.glm.2)$dispersion, n) #rep(MASS::gamma.shape(start.temp.glm.2)$alpha, n)
  alpha3.current <- rep(1/summary(start.temp.glm.3)$dispersion, n) #rep(MASS::gamma.shape(start.temp.glm.3)$alpha, n)
  beta.current.1 <- 1/summary(start.temp.glm.1)$dispersion / predict(start.temp.glm.1, type="response")#, dispersion = 1/MASS::gamma.shape(start.temp.glm.1) )
  beta.current.2 <- 1/summary(start.temp.glm.2)$dispersion / predict(start.temp.glm.2, type="response")#, dispersion = 1/MASS::gamma.shape(start.temp.glm.2) )
  beta.current.3 <- 1/summary(start.temp.glm.3)$dispersion / predict(start.temp.glm.3, type="response")#, dispersion = 1/MASS::gamma.shape(start.temp.glm.3) )
  beta.current   <- rowMeans(cbind(beta.current.1, beta.current.2, beta.current.3))
  #--------------------------

  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)
  loglike.store   <- rep(-.Machine$double.xmax, 3)
  if (is.null(l3)) stop("l3 is NULL, Double Gamma Regression")
  
  coef1.trace  <- coef1.current
  coef2.trace  <- coef2.current
  coef3.trace  <- coef3.current
  alpha1.trace <- alpha1.current
  alpha2.trace <- alpha2.current
  alpha3.trace <- alpha3.current
  beta.trace   <- beta.current
  j            <- 1
  
  while ( (loglike.diff > tol) && (j <= maxit) ) {
    
    #-------
    #E step
    #-------
    Expected.x3    <- rep(0,n)
    Expected.x1    <- rep(0,n)
    Expected.x2    <- rep(0,n)
    Expected.logx1 <- rep(0,n)
    Expected.logx2 <- rep(0,n)
    Expected.logx3 <- rep(0,n)
    den            <- rep(0,n)
    for (i in seq_len(n)){
      expected.latent.res <- expected.latent(c(y1[i],y2[i]),
                                             alpha=c(alpha1.current[i],alpha2.current[i],alpha3.current[i]),
                                             beta=beta.current[i])
      Expected.x3[i]      <- expected.latent.res$Expected.s
      Expected.x1[i]      <- y1[i]-expected.latent.res$Expected.s
      Expected.x2[i]      <- y2[i]-expected.latent.res$Expected.s
      Expected.logx3[i]   <- expected.latent.res$Expected.logs
      Expected.logx1[i]   <- expected.latent.res$Expected.logy1s
      Expected.logx2[i]   <- expected.latent.res$Expected.logy2s
      if (Expected.x1[i] <=0 || Expected.x2[i] <=0){
          Expected.x3[i]  <- min(y1[i], y2[i])/2
          Expected.x1[i]  <- y1[i]-Expected.x3[i]
          Expected.x2[i]  <- y2[i]-Expected.x3[i]
          warning("negative expected x1 or x2 at ", i, " obs at ", j, " iteration", "\n")
      }
      den[i]              <- expected.latent.res$denominator
      #if (expected.latent.res$numerator.integration.message != "OK"){print(expected.latent.res$numerator.integration.message)}
      #if (expected.latent.res$denominator.integration.message != "OK"){print(expected.latent.res$denominator.integration.message)}
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
    
    Q.function <- function(coef, Model.Matrix, Expected.log.value){
      q.res <- sum(exp(coef%*%t(Model.Matrix)) * log(beta.current)) - 
        sum(lgamma(exp(coef%*%t(Model.Matrix)))) +
        sum(exp(coef%*%t(Model.Matrix))*Expected.log.value)
      return(q.res)
    }
    
    # method 1 # using optim instead of glm function
    # data1$x1  <- alpha1.current
    # GLM.loglike <- function(theta, Model.Matrix, data){
    #   E <- exp(theta[-1]%*%t(Model.Matrix))
    #   V <- exp(theta[1])*E^2
    #   glm.loglik <- sum(dgamma(data$x1, shape=1/exp(theta[1]), rate=1/(E*exp(theta[1])), log=TRUE))
    #   return(glm.loglik)
    # }
    # coef1.new.temp <- optim(par=as.vector(c(log(dispersion1.current), coef1.current)), 
    #                         fn=GLM.loglike, Model.Matrix=Model.Matrix.1, data=data1, 
    #                         gr=NULL, hessian = TRUE,
    #                         control = list(fnscale = -1))
    # #trace(GLM.loglike, tracer=quote(print(theta)))
    # dispersion1.new <- exp(coef1.new.temp$par[1])
    # method 2 # 
    # this gives a good starting value for the optim function later
    m1             <- glm(alpha1.current~Model.Matrix.1[,-1], family=Gamma(link="log"))
    coef1.new.temp <- as.vector(m1$coefficients)
    coef1.new      <- optim(par=coef1.new.temp, 
                            fn=Q.function, 
                            Model.Matrix=Model.Matrix.1, Expected.log.value=Expected.logx1,
                            gr=NULL, hessian=TRUE, control=list(fnscale=-1))$par
    alpha1.new     <- as.vector(exp(coef1.new%*%t(Model.Matrix.1)))  
    
    m2             <- glm(alpha2.current~Model.Matrix.2[,-1], family=Gamma(link="log"))
    coef2.new.temp <- as.vector(m2$coefficients)
    coef2.new      <- optim(par=coef2.new.temp, 
                            fn=Q.function, 
                            Model.Matrix=Model.Matrix.2, Expected.log.value=Expected.logx2,
                            gr=NULL, hessian=TRUE, control=list(fnscale=-1))$par
    alpha2.new     <- as.vector(exp(coef2.new%*%t(Model.Matrix.2)))
    
    m3             <- glm(alpha3.current~Model.Matrix.3[,-1], family=Gamma(link="log"))
    coef3.new.temp <- as.vector(m3$coefficients)
    coef3.new      <- optim(par=coef3.new.temp, 
                            fn=Q.function, 
                            Model.Matrix=Model.Matrix.3, Expected.log.value=Expected.logx3,
                            gr=NULL, hessian=TRUE, control=list(fnscale=-1))$par
    alpha3.new     <- as.vector(exp(coef3.new%*%t(Model.Matrix.3)))
    
    beta.new <- (alpha1.new + alpha2.new + alpha3.new) / (Expected.x1 + Expected.x2 + Expected.x3)
    
    if (verbose){
      cat(c("alpha1:", mean(alpha1.new), "\n"))
      cat(c("alpha2:", mean(alpha2.new), "\n"))
      cat(c("alpha3:", mean(alpha3.new), "\n"))
      cat(c("beta:", mean(beta.new), "\n"))
    }
    
    coef1.current  <- coef1.new
    coef2.current  <- coef2.new
    coef3.current  <- coef3.new
    
    loglike.current<- loglike.new
    
    alpha1.current <- alpha1.new 
    alpha2.current <- alpha2.new 
    alpha3.current <- alpha3.new 
    beta.current   <- beta.new 
    
    alpha1.trace   <- cbind(alpha1.trace, alpha1.new)
    alpha2.trace   <- cbind(alpha2.trace, alpha2.new)
    alpha3.trace   <- cbind(alpha3.trace, alpha3.new)
    beta.trace     <- cbind(beta.trace, beta.new)
    coef1.trace    <- cbind(coef1.trace, coef1.new)
    coef2.trace    <- cbind(coef2.trace, coef2.new)
    coef3.trace    <- cbind(coef3.trace, coef3.new)
    
    j              <- j+1
  
  }
  
  if (!all(loglike[1:(j-1)] == cummax(loglike[1:(j-1)]))){
    cat("Loglikelihood is not monotonically increasing!", "\n")
  }
  
  
  noparams<- m1$rank + m2$rank + m3$rank + n
  AIC     <- -2*loglike[j-1] + noparams * 2
  BIC     <- -2*loglike[j-1] + noparams * log(n)
  
  # fitted values
  y1.fitted <- (alpha1.current + alpha3.current) / beta.current
  y2.fitted <- (alpha2.current + alpha3.current) / beta.current
  y1.residual <- y1 - y1.fitted
  y2.residual <- y2 - y2.fitted
  
  # formula
  l1.formula <- formula(paste(namey1, as.character(l1.new)[3], sep="~"))
  l2.formula <- formula(paste(namey2, as.character(l2.new)[3], sep="~"))
  l3.formula <- formula(paste("", as.character(l2.new)[3], sep="~"))
  
  result<-list(coefficients=list(coef1.new, coef2.new, coef3.new),
               alpha1      =alpha1.current,
               alpha2      =alpha2.current,
               alpha3      =alpha3.current,
               beta        =beta.current,
               
               alpha1.trace=alpha1.trace,
               alpha2.trace=alpha2.trace,
               alpha3.trace=alpha3.trace,
               beta.trace  =beta.trace,
               coef.trace  =coef3.trace,
               
               fitted      =data.frame(y1=y1.fitted,
                                       y2=y2.fitted),
               loglikelihood=loglike[j-1],
               llseq       = loglike[1:(j-1)],
               
               parameters.number=noparams,
               AIC         = AIC,
               BIC         = BIC,
               
               y           = cbind(y1, y2),
               Model.Matrix= list(Model.Matrix.1, Model.Matrix.2, Model.Matrix.3),
               formula     = list(l1.formula, l2.formula, l3.formula),
               
               #call=tempcall,
               iterations=j-1)
  options(warn=0)
  class(result)<-c('BGR', 'glm')
  return(result)
}

#
try.aviva <- BGR(data   = nnaviva, # data = naviva,
           l1     = COST_AD ~ PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor + PH_PENPTS + MD_LIC_CAT_A,
           l2     = COST_TPD ~ PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor+ PH_PENPTS + MD_LIC_CAT_A,
           l3     = ~ PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor+ PH_PENPTS + MD_LIC_CAT_A,
           expo   = NULL,
           maxit  = 100,
           tol    = 1e-4,
           Aitken = FALSE,
           verbose= TRUE)
all(try.aviva$llseq == cummax(try.aviva$llseq))
plot(try.aviva$llseq)
summary(try.aviva$alpha1) 
summary(try.aviva$alpha2)
summary(try.aviva$alpha3)
summary(try.aviva$beta)
boxplot(try.aviva$alpha1, ylim=c(0, 15), main="alpha1") 
boxplot(try.aviva$alpha2, ylim=c(0, 15), main="alpha2")
boxplot(try.aviva$alpha3, ylim=c(0, 15), main="alpha3")
boxplot(try.aviva$beta, ylim=c(0, 15), main="beta")

plot(nnaviva$COST_AD, nnaviva$COST_TPD, main="by mvGamma")
points(try.aviva$fitted, col=2, cex=1, pch=20)

uni.mod1 <- glm(COST_AD ~ PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor + PH_PENPTS + MD_LIC_CAT_A,
                family=Gamma(link="log"), data=nnaviva)
mod1.pred <- predict(uni.mod1, type="response", dispersion=1/MASS::gamma.shape(uni.mod1))
uni.mod2 <- glm(COST_TPD ~ PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor + PH_PENPTS + MD_LIC_CAT_A,
                family=Gamma(link="log"), data=nnaviva)
mod2.pred <- predict(uni.mod2, type="response", dispersion=1/MASS::gamma.shape(uni.mod2))   #uni.mod1$fitted.values
plot(nnaviva$COST_AD, nnaviva$COST_TPD, main="by independent models")
points(mod1.pred, mod2.pred, col=4, cex=1, pch=20)

round(ComparisonMetrics(nnaviva$COST_AD, try.aviva$fitted[,1], expo = rep(1, 2074)),4)
round(ComparisonMetrics(nnaviva$COST_AD, mod1.pred, expo = rep(1, 2074)), 4)
summary(nnaviva$COST_AD); summary(try.aviva$fitted[,1]); summary(mod1.pred)

round(ComparisonMetrics(nnaviva$COST_TPD, try.aviva$fitted[,2], expo = rep(1, 2074)),4)
round(ComparisonMetrics(nnaviva$COST_TPD, mod2.pred, expo = rep(1, 2074)), 4)
summary(nnaviva$COST_TPD); summary(try.aviva$fitted[,2]); summary(mod2.pred)


#----------
# artifitial data

common <- rgamma(500, shape = 4, rate=.2)
response.var1 <- rgamma(500, shape = 5, rate=.2) + common
response.var2 <- rgamma(500, shape = 3, rate=.1) + common
plot(response.var1, response.var2)
cor(response.var1, response.var2)
x1 <- (response.var1)/2+(response.var2)^(1/3)+rgamma(500, shape=3, rate=.4) 
x2 <- (response.var1)/4+sqrt(response.var2)+rgamma(500, shape=4, rate=.6) 
plot(x1, x2)

arti.data <- as.data.frame(cbind(response.var1, response.var2, x1, x2))
names(arti.data)
names(arti.data)[3] <- "predictor.x1"
names(arti.data)[4] <- "predictor.x2"

summary(g1 <- glm(response.var1 ~ predictor.x1+predictor.x2, data=arti.data, family=Gamma(link="log")))
summary(g2 <- glm(response.var2 ~ predictor.x1+predictor.x2, data=arti.data, family=Gamma(link="log")))
g1.pred <- predict(g1, type="response", dispersion=1/MASS::gamma.shape(g1))
g2.pred <- predict(g2, type="response", dispersion=1/MASS::gamma.shape(g2))
plot(response.var1, g1.pred)
plot(response.var2, g2.pred)
plot(g1.pred, g2.pred)


###########
#arti.data.save <- arti.data
###########

try <- BGR(data   = arti.data,
           l1     = response.var1 ~ predictor.x1 + predictor.x2,
           l2     = response.var2 ~ predictor.x1 + predictor.x2,
           l3     = ~ predictor.x1 + predictor.x2,
           expo   = NULL,
           maxit  = 400,
           tol    = 1e-5,
           Aitken = FALSE,
           verbose= TRUE)
plot(try$llseq)
all(try$llseq == cummax(try$llseq))
summary(try$alpha1)
summary(try$alpha2)
summary(try$alpha3)
summary(try$beta)
summary(try$fitted)
plot(arti.data$response.var1, arti.data$response.var2, main="by mvGamma")
points(try$fitted, col=2, cex=1, pch=20)
plot(arti.data$response.var1, arti.data$response.var2, main="by independent models")
points(g1.pred, g2.pred, col=4, cex=1, pch=20)


round(ComparisonMetrics(arti.data$response.var1, try$fitted[,1], expo = rep(1, 500)),4)
round(ComparisonMetrics(arti.data$response.var1, g1.pred, expo = rep(1, 500)), 4)
summary(arti.data$response.var1); summary(try$fitted[,1]); summary(g1.pred)

round(ComparisonMetrics(arti.data$response.var2, try$fitted[,2], expo = rep(1, 500)), 4)
round(ComparisonMetrics(arti.data$response.var2, g2.pred, expo = rep(1, 500)), 4)
summary(arti.data$response.var2); summary(try$fitted[,2]); summary(g2.pred)


plot(colMeans( try$alpha1.trace))




aviva1 <- BGR(data   = naviva, # data = naviva,
              l1     = COST_AD ~ PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor + PH_PENPTS + MD_LIC_CAT_A,
              l2     = COST_TPD ~ PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor+ PH_PENPTS + MD_LIC_CAT_A,
              l3     = ~ PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor+ PH_PENPTS + MD_LIC_CAT_A,
              expo   = NULL,
              maxit  = 200,
              tol    = 1e-5,
              Aitken = FALSE,
              verbose= TRUE)
plot(aviva1$llseq)
aviva1$BIC
aviva1$llseq[28]

train.index <- sample(c(1:n), 1450, replace = FALSE)
naviva.train <- naviva[train.index,]
naviva.test  <- naviva[setdiff(1:n, train.index),]

aviva1.train <- BGR(data   = naviva.train, # data = naviva,
                    l1     = COST_AD ~ PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor + PH_PENPTS + MD_LIC_CAT_A,
                    l2     = COST_TPD ~ PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor+ PH_PENPTS + MD_LIC_CAT_A,
                    l3     = ~ PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor+ PH_PENPTS + MD_LIC_CAT_A,
                    expo   = NULL,
                    maxit  = 200,
                    tol    = 1e-5,
                    Aitken = FALSE,
                    verbose= TRUE)

aviva1.train.copy <- BGR(data   = naviva.train, # data = naviva,
                         l1     = COST_AD ~ PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor + PH_PENPTS + MD_LIC_CAT_A,
                         l2     = COST_TPD ~ PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor+ PH_PENPTS + MD_LIC_CAT_A,
                         l3     = ~ PH_GEND + VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK + EXCESS_TOT.factor+ PH_PENPTS + MD_LIC_CAT_A,
                         expo   = NULL,
                         maxit  = 200,
                         tol    = 1e-5,
                         Aitken = FALSE,
                         verbose= TRUE)
