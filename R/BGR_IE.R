

BGR.IE <- function(data,
                   response,
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
  data$weight.var <- weight.var

  data1 <- data;
  data1 <- data1[,names(data1)!=namey1];
  data1 <- data1[,names(data1)!=namey2]
  datap <- data1

  if (!is.null(l4)){
    if (as.character(l4[2])=="."){l4.new <- formula( paste( "be", paste( names(data1), "", collapse = "+", sep = ""), sep="~"))}
    else {l4.new <- formula(paste("be", as.character(l4)[2], sep="~"))}
  } else stop("l4 cannot be NULL, must be supplied")

  #-----------------------------
  # starting values
  x3.start <- rep(0, n) # set.seed(10)
  for (i in seq_len(n)){
    x3.start[i] <- runif(1, min=0, max=min(y1[i], y2[i]))
  }
  x1.start         <- y1 - x3.start
  x2.start         <- y2 - x3.start

  alpha1.current   <- MASS::fitdistr(x1.start, "gamma")$estimate[1]
  alpha2.current   <- MASS::fitdistr(x1.start, "gamma")$estimate[1]
  alpha3.current   <- MASS::fitdistr(x1.start, "gamma")$estimate[1]
  beta.current     <- rowMeans(cbind(alpha1.current/x1.start, alpha2.current/x2.start, alpha3.current/x3.start))

  data.start       <- data1
  data.start$be    <- beta.current
  start.temp.glm   <- glm(l4.new, data=data.start, family=Gamma(link="log"))
  coef4.current    <- start.temp.glm$coefficient
  Model.Matrix.4   <- model.matrix(start.temp.glm)

  #--------------------------

  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.eps)
  loglike         <- rep(0,maxit)
  loglike.store   <- rep(-.Machine$double.xmax, 3)

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
                                             alpha=c(alpha1.current,alpha2.current,alpha3.current),
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
    Q.b.function   <- function(coef){
      q.b.res <- sum((alpha1.current + alpha2.current + alpha3.current) * log(exp(coef %*% t(Model.Matrix.4)))) -
        sum(exp(coef %*% t(Model.Matrix.4)) * (Expected.x1 + Expected.x2 + Expected.x3))
      return(q.b.res)
    }
    beta.temp      <- (alpha1.current+alpha2.current+alpha3.current) / (Expected.x1+Expected.x2+Expected.x3)
    m4             <- glm(beta.temp~Model.Matrix.4[,-1], family=Gamma(link="log"))
    coef4.new      <- optim(par     = as.vector(m4$coefficients),
                            fn      = Q.b.function,
                            gr      = NULL,
                            hessian = TRUE,
                            control = list(fnscale=-1))$par
    beta.new       <- as.vector(exp(coef4.new%*%t(Model.Matrix.4)))

    alpha1.rootfun <- function(alpha1.var){
      sum(log(beta.new)) - digamma(alpha1.var)*n + sum(Expected.logx1)
    }
    alpha1.new <- uniroot(alpha1.rootfun,
                          lower=sqrt(.Machine$double.eps),
                          upper=100000,
                          tol = sqrt(.Machine$double.xmin),
                          check.conv = TRUE)$root

    alpha2.rootfun <- function(alpha2.var){
      sum(log(beta.new)) - digamma(alpha2.var)*n + sum(Expected.logx2)
    }
    alpha2.new <- uniroot(alpha2.rootfun,
                          lower=sqrt(.Machine$double.eps),
                          upper=100000,
                          tol = sqrt(.Machine$double.xmin),
                          check.conv = TRUE)$root
    alpha3.rootfun <- function(alpha3.var){
      sum(log(beta.new)) - digamma(alpha3.var)*n + sum(Expected.logx3)
    }
    alpha3.new <- uniroot(alpha3.rootfun,
                          lower=sqrt(.Machine$double.eps),
                          upper=100000,
                          tol = sqrt(.Machine$double.xmin),
                          check.conv = TRUE)$root

    if (verbose){
      cat(c("alpha1:", round(alpha1.new,4), "\n"))
      cat(c("alpha2:", round(alpha2.new,4), "\n"))
      cat(c("alpha3:", round(alpha3.new,4), "\n"))
      cat(c("beta:",   round(summary(beta.new),4),  "\n"))
    }

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
    coef4.trace    <- cbind(coef4.trace, coef4.new)
    j              <- j+1
  }

  if (!all(loglike[1:(j-1)] == cummax(loglike[1:(j-1)]))){
    cat("Loglikelihood is not monotonically increasing!", "\n")
  }


  noparams<- m4$rank + 3
  AIC     <- -2*loglike[j-1] + noparams * 2
  BIC     <- -2*loglike[j-1] + noparams * log(n)

  # fitted values
  y1.fitted <- (alpha1.current + alpha3.current) / beta.current
  y2.fitted <- (alpha2.current + alpha3.current) / beta.current
  y1.residual <- y1 - y1.fitted
  y2.residual <- y2 - y2.fitted

  # formula
  l4.formula <- formula(paste("", as.character(l4.new)[3], sep="~"))

  result<-list(coefficients=coef4.new,
               alpha1      =alpha1.current,
               alpha2      =alpha2.current,
               alpha3      =alpha3.current,
               beta        =beta.current,

               alpha1.trace=alpha1.trace,
               alpha2.trace=alpha2.trace,
               alpha3.trace=alpha3.trace,
               beta.trace  =beta.trace,
               coef.trace  =coef4.trace,

               fitted      =data.frame(y1=y1.fitted,
                                       y2=y2.fitted),
               loglikelihood=loglike[j-1],
               llseq       = loglike[1:(j-1)],

               parameters.number=noparams,
               AIC         = AIC,
               BIC         = BIC,

               y           = cbind(y1, y2),
               Model.Matrix= list(Model.Matrix.4),
               formula     = list(l4.formula),

               call        = tempcall,
               iterations  = j-1)
  options(warn=0)
  class(result)<-c('BGR.IE', 'glm')
  return(result)
}

