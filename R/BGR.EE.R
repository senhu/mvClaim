
BGR.EE <- function(data,
                   response,
                   l1,
                   l2,
                   l3,
                   l4,
                   expo   = NULL,
                   maxit  = 200,
                   tol    = 1e-4,
                   Aitken = FALSE,
                   verbose= TRUE){
  #------------------------------
  options(warn=-1)
  tempcall <- as.call( c(expression(BGR), list(data    = substitute(data),
                                               response= response,
                                               l1      = l1,
                                               l2      = l2,
                                               l3      = l3,
                                               l4      = l4,
                                               expo    = expo,
                                               maxit   = maxit,
                                               tol     = tol,
                                               Aitken  = Aitken,
                                               verbose = verbose) ) )
  #------------------------------
  namey1 <- response[1]
  namey2 <- response[2]
  y1     <- data[,names(data)==namey1]
  y2     <- data[,names(data)==namey2]
  n      <- length(y1)

  if (!is.null(expo)) {weight.var<-data[,names(data)==expo]} else {weight.var<-rep(1, n)}
  data$weight.var <- weight.var # there should be better ways of doing this

  data1 <- data;
  data1 <- data1[, names(data1)!=namey1];
  data1 <- data1[, names(data1)!=namey2]
  datap <- data1

  if (!is.null(l1)){
    if (as.character(l1[2])=="."){
      l1.new <- formula( paste( "x1", paste( names(data1),'',collapse='+',sep='' ), sep='~'))
    } else {l1.new <- formula(paste("x1", as.character(l1)[2], sep="~"))}
  } else stop("l1 cannot be NULL, must be supplied")
  if (!is.null(l2)){
    if (as.character(l2[2])=="."){l2.new <- formula( paste( "x2", paste( names(data1),'',collapse='+',sep='' ), sep='~'))}
    else {l2.new <- formula(paste("x2", as.character(l2)[2], sep="~"))}
  } else stop("l2 cannot be NULL, must be supplied")
  if (!is.null(l3)){
    if (as.character(l3[2])=="."){l3.new <- formula( paste( "x3", paste( names(data1), "", collapse = "+", sep = ""), sep="~"))}
    else {l3.new <- formula(paste("x3", as.character(l3)[2], sep="~"))}
  } else stop("l3 cannot be NULL, must be supplied")
  if (!is.null(l4)){
    if (as.character(l4[2])=="."){l4.new <- formula( paste( "be", paste( names(data1), "", collapse = "+", sep = ""), sep="~"))}
    else {l4.new <- formula(paste("be", as.character(l4)[2], sep="~"))}
  } else stop("l4 cannot be NULL, must be supplied")

  #-----------------------------
  # starting values
  x3.start <- rep(0, n)
  #set.seed(10)
  for (i in seq_len(n)){
    x3.start[i] <- runif(1, min=0, max=min(y1[i], y2[i]))
  }
  x1.start         <- y1 - x3.start
  x2.start         <- y2 - x3.start
  data.start       <- data1
  data.start$x1    <- x1.start
  data.start$x2    <- x2.start
  data.start$x3    <- x3.start
  start.temp.glm.1 <- glm(l1.new, data=data.start, family=Gamma(link="log"))
  start.temp.glm.2 <- glm(l2.new, data=data.start, family=Gamma(link="log"))
  start.temp.glm.3 <- glm(l3.new, data=data.start, family=Gamma(link="log"))

  Model.Matrix.1   <- model.matrix(start.temp.glm.1)
  Model.Matrix.2   <- model.matrix(start.temp.glm.2)
  Model.Matrix.3   <- model.matrix(start.temp.glm.3)
  coef1.current    <- start.temp.glm.1$coefficient
  coef2.current    <- start.temp.glm.2$coefficient
  coef3.current    <- start.temp.glm.3$coefficient

  alpha1.current   <- rep(1/summary(start.temp.glm.1)$dispersion, n) #rep(MASS::gamma.shape(start.temp.glm.1)$alpha, n)
  alpha2.current   <- rep(1/summary(start.temp.glm.2)$dispersion, n) #rep(MASS::gamma.shape(start.temp.glm.2)$alpha, n)
  alpha3.current   <- rep(1/summary(start.temp.glm.3)$dispersion, n) #rep(MASS::gamma.shape(start.temp.glm.3)$alpha, n)

  beta.current.1   <- 1/summary(start.temp.glm.1)$dispersion / predict(start.temp.glm.1, type="response")#, dispersion = 1/MASS::gamma.shape(start.temp.glm.1) )
  beta.current.2   <- 1/summary(start.temp.glm.2)$dispersion / predict(start.temp.glm.2, type="response")#, dispersion = 1/MASS::gamma.shape(start.temp.glm.2) )
  beta.current.3   <- 1/summary(start.temp.glm.3)$dispersion / predict(start.temp.glm.3, type="response")#, dispersion = 1/MASS::gamma.shape(start.temp.glm.3) )
  beta.current     <- rowMeans(cbind(beta.current.1, beta.current.2, beta.current.3))

  data.start$be    <- beta.current
  start.temp.glm.4 <- glm(l4.new, data=data.start, family=Gamma(link="log"))
  coef4.current    <- start.temp.glm.4$coefficient
  Model.Matrix.4   <- model.matrix(start.temp.glm.4)

  # p1 <- dim(Model.Matrix.1)[2]
  # p2 <- dim(Model.Matrix.2)[2]
  # p3 <- dim(Model.Matrix.3)[2]
  # p4 <- dim(Model.Matrix.4)[2]
  #--------------------------
  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)
  loglike.store   <- rep(-.Machine$double.xmax, 3)
  #if (is.null(l3)) stop("l3 is NULL, Double Gamma Regression")

  coef1.trace  <- coef1.current
  coef2.trace  <- coef2.current
  coef3.trace  <- coef3.current
  coef4.trace  <- coef4.current
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
                                             beta =beta.current[i])
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
      cat(c("iter:", j, "; current loglike:", round(loglike.new,4), "; loglike change:", round(loglike.diff,6), "\n"))
    }

    #--------
    # M step
    #--------

    beta.temp      <- (alpha1.current + alpha2.current + alpha3.current) / (Expected.x1 + Expected.x2 + Expected.x3)
    m4             <- glm(beta.temp ~ Model.Matrix.4[,-1], family=Gamma(link="log"))
    Q.b.function   <- function(coef){
      q.b.res <- sum((alpha1.current + alpha2.current + alpha3.current) * log(exp(coef %*% t(Model.Matrix.4)))) -
        sum(exp(coef %*% t(Model.Matrix.4)) * (Expected.x1 + Expected.x2 + Expected.x3))
      return(q.b.res)
    }
    coef4.optim    <- optim(par     = as.vector(m4$coefficients),
                            fn      = Q.b.function,
                            gr      = NULL,
                            hessian = TRUE,
                            control = list(fnscale=-1, maxit=1e+5))
    if (coef4.optim$convergence!=0) cat("optim coef. beta not converged")
    coef4.new      <- coef4.optim$par
    beta.new       <- as.vector(exp(coef4.new%*%t(Model.Matrix.4)))

    m1             <- glm(alpha1.current~Model.Matrix.1[,-1], family=Gamma(link="log"))
    Q.function.a1  <- function(coef){
      q.res <- sum( exp(coef %*% t(Model.Matrix.1)) * log(beta.new) ) -
               sum( lgamma(exp(coef%*%t(Model.Matrix.1))) ) +
               sum( exp(coef %*% t(Model.Matrix.1)) * Expected.logx1 )
      return(q.res)
    }
    coef1.optim    <- optim(par     = as.vector(m1$coefficients),
                            fn      = Q.function.a1,
                            gr      = NULL,
                            hessian = TRUE,
                            control = list(fnscale=-1, maxit=1e+5))
    if (coef1.optim$convergence!=0) cat("optim coef a1 not converged", "\n")
    coef1.new      <- coef1.optim$par
    alpha1.new     <- as.vector(exp(coef1.new%*%t(Model.Matrix.1)))

    Q.function.a2  <- function(coef){
      q.res <- sum( exp(coef %*% t(Model.Matrix.2)) * log(beta.new) ) -
               sum( lgamma(exp(coef%*%t(Model.Matrix.2))) ) +
               sum( exp(coef %*% t(Model.Matrix.2)) * Expected.logx2 )
      return(q.res)
    }
    m2             <- glm(alpha2.current~Model.Matrix.2[,-1], family=Gamma(link="log"))
    coef2.optim    <- optim(par     = as.vector(m2$coefficients),
                            fn      = Q.function.a2,
                            gr      = NULL,
                            hessian = TRUE,
                            control = list(fnscale=-1, maxit=1e+5))
    if (coef2.optim$convergence) cat("optim coef a2 not converged", "\n")
    coef2.new      <- coef2.optim$par
    alpha2.new     <- as.vector(exp(coef2.new%*%t(Model.Matrix.2)))

    Q.function.a3 <- function(coef){
      q.res <- sum( exp(coef %*% t(Model.Matrix.3)) * log(beta.new) ) -
               sum( lgamma(exp(coef%*%t(Model.Matrix.3))) ) +
               sum( exp(coef %*% t(Model.Matrix.3)) * Expected.logx3 )
      return(q.res)
    }
    m3             <- glm(alpha3.current~Model.Matrix.3[,-1], family=Gamma(link="log"))
    coef3.optim    <- optim(par     = as.vector(m3$coefficients),
                            fn      = Q.function.a3,
                            gr      = NULL,
                            hessian = TRUE,
                            control = list(fnscale=-1, maxit=1e+5))
    if (coef3.optim$convergence!=0) cat("optim coef a3 not converged", "\n")
    coef3.new <- coef3.optim$par
    alpha3.new     <- as.vector(exp(coef3.new%*%t(Model.Matrix.3)))

    if (verbose){
      cat(c("alpha1:", round(summary(alpha1.new),4), "\n"))
      cat(c("alpha2:", round(summary(alpha2.new),4), "\n"))
      cat(c("alpha3:", round(summary(alpha3.new),4), "\n"))
      cat(c("beta:",   round(summary(beta.new),4),  "\n"))
    }

    coef1.current  <- coef1.new
    coef2.current  <- coef2.new
    coef3.current  <- coef3.new
    coef4.current  <- coef4.new
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
    coef4.trace    <- cbind(coef4.trace, coef4.new)
    j              <- j+1
  }

  if (!all(loglike[1:(j-1)] == cummax(loglike[1:(j-1)]))){
    cat("Loglikelihood is not monotonically increasing!", "\n")
  }

  noparams<- m1$rank + m2$rank + m3$rank + m3$rank
  AIC     <- -2*loglike[j-1] + noparams * 2
  BIC     <- -2*loglike[j-1] + noparams * log(n)

  # fitted values
  y1.fitted <- (alpha1.current + alpha3.current) / beta.current
  y2.fitted <- (alpha2.current + alpha3.current) / beta.current
  y1.residual <- y1 - y1.fitted
  y2.residual <- y2 - y2.fitted

  # formula
  l1.formula <- formula(paste("", as.character(l1.new)[3], sep="~"))
  l2.formula <- formula(paste("", as.character(l2.new)[3], sep="~"))
  l3.formula <- formula(paste("", as.character(l3.new)[3], sep="~"))
  l4.formula <- formula(paste("", as.character(l4.new)[3], sep="~"))

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
               formula     = list(l1.formula, l2.formula, l3.formula, l4.formula),

               call        = tempcall,
               iterations  = j-1)
  options(warn=0)
  class(result)<-c('BGR.EE', 'glm')
  return(result)
}



