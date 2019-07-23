#' Bivariate Gamma Distribution Estimation
#'
#' This function allows you to estimate bivariate gamma distribution given a data set using EM algorithm
#'
#' @param data A numeric vector, matrix or data frame of observations. Categorical variables are not allowed. If a matrix or data frame, rows correspond to observations and columns correspond to variables.
#' @param maxit An integer limits on the number of EM iterations. The default is 300.
#' @param tol A value giving relative convergence tolerance for the log-likelihood. The default is 1e-6.
#' @param start Starting value of the EM algorithm
#' @param verbose logical; controls whether summary result of each EM iteration is displayed during the fitting procedure. Default is TRUE.
#'
#' @return An object of class \code{BGR} providing the estimation results.
#'   The details of the output components are:
#'   \item{estimate}{The estimated parametr values.}
#'   \item{loglike}{The final estimated maximum log-likelihood value.}
#'   \item{ll}{The sequence of log-likelihood values in the EM algorithm fitting process.}
#'   \item{df}{The number of estimated parameters.}
#'   \item{AIC}{AIC values.}
#'   \item{BIC}{BIC values.}
#'   \item{iter}{Total iteration numbers.}
#'   \item{n}{The number of observations in the data.}
#'   \item{call}{The matched call.}
#'
#' @example
#' ndat <- rbivgamma(500, alpha = c(1,2,0.5), beta=0.1)
#' plot(ndat, xlab="", ylab="")
#' mod <- BGE(data = ndat,
#'            maxit = 100, tol = 1e-5,
#'            verbose = TRUE)
#' mod
#'


BGE <- function(data,
                maxit = 200,
                tol = 1e-5,
                start = NULL,
                verbose = TRUE)
{
  tempcall<-as.call(c(expression(BGE), list(data    = substitute(data),
                                            maxit   = maxit,
                                            tol     = tol,
                                            start   = start,
                                            verbose = verbose)) )

  n<-dim(data)[1]
  p<-dim(data)[2]
  if (p!=2){stop("The data set is not bivariate.")}
  y1 <- data[,1]
  y2 <- data[,2]

  # starting value
  if (is.null(start)){
    t3 <- apply(data, 1, min)/2
    t.y1 <- y1-t3
    t.y2 <- y2-t3
    fitd.1 <- MASS::fitdistr(t.y1/mean(t.y1), "gamma")
    fitd.2 <- MASS::fitdistr(t.y2/mean(t.y2), "gamma")
    fitd.3 <- MASS::fitdistr(t3/mean(t3), "gamma")
    beta.current <- mean(c(fitd.1$estimate[2]/mean(t.y1),
                           fitd.2$estimate[2]/mean(t.y2),
                           fitd.3$estimate[2]/mean(t3)))
    alpha1.current <- fitd.1$estimate[1]
    alpha2.current <- fitd.2$estimate[1]
    alpha3.current <- fitd.3$estimate[1]
  } else {
    if (length(start) != 4){stop("wrong starting values provided")}
    alpha1.current <- start[1]
    alpha2.current <- start[2]
    alpha3.current <- start[3]
    beta.current   <- start[4]
  }

  j<-1
  loglike.diff    <- 1000
  loglike         <- rep(0,maxit)
  loglike.store   <- rep(-sqrt(.Machine$double.xmax), 3)

  loglike.current <- BGE.data.loglikelihood(data,
                                            alpha=c(alpha1.current, alpha2.current, alpha3.current),
                                            beta=beta.current)
  if (verbose){
    cat(c("starting loglikelihood:", round(loglike.current,5), "\n"))
  }
  while ( (loglike.diff >= tol) && (j <= maxit) ) {
    #-------
    #E step
    #-------
    Estep.res <- BGE.Estep.function(data,
                                    alpha=c(alpha1.current, alpha2.current, alpha3.current),
                                    beta=beta.current)
    Expected.x3 <- Estep.res[[1]]
    Expected.x1 <- (data[,1]-Expected.x3)
    Expected.x2 <- (data[,2]-Expected.x3)
    Expected.logx3 <- Estep.res[[2]]
    Expected.logx1 <- Estep.res[[3]]
    Expected.logx2 <- Estep.res[[4]]

    #-------
    #M step
    #-------
    beta.new <- ((alpha1.current+alpha2.current+alpha3.current)*n) / sum(c(Expected.x1,Expected.x2,Expected.x3))

    alpha1.rootfun <- function(alpha1.var){
      log(beta.new)*n - digamma(alpha1.var)*n + sum(Expected.logx1)
    }
    alpha1.new <- stats::uniroot(alpha1.rootfun,
                          lower=sqrt(.Machine$double.eps),
                          upper=100000,
                          tol = sqrt(.Machine$double.xmin),
                          check.conv = TRUE)$root
    alpha2.rootfun <- function(alpha2.var){
      log(beta.new)*n - digamma(alpha2.var)*n + sum(Expected.logx2)
    }
    alpha2.new <- stats::uniroot(alpha2.rootfun,
                          lower=sqrt(.Machine$double.eps),
                          upper=100000,
                          tol = sqrt(.Machine$double.xmin),
                          check.conv = TRUE)$root
    alpha3.rootfun <- function(alpha3.var){
      log(beta.new)*n - digamma(alpha3.var)*n + sum(Expected.logx3)
    }
    alpha3.new <- stats::uniroot(alpha3.rootfun,
                          lower=sqrt(.Machine$double.eps),
                          upper=100000,
                          tol=sqrt(.Machine$double.xmin),
                          check.conv = TRUE)$root

    loglike.new <- BGE.data.loglikelihood(data,
                                          alpha=c(alpha1.new, alpha2.new, alpha3.new),
                                          beta=beta.new)

    if (verbose){
      cat(c("alpha1:", round(alpha1.new,4),"alpha2:", round(alpha2.new,4),"alpha3:", round(alpha3.new,4), "\n"))
      cat(c("beta:", round(beta.new,4), "\n"))
    }

    loglike.store <- c(loglike.store[-1], loglike.new)
    loglike.diff <- abs(loglike.new - loglike.current) / (1+abs(loglike.new))

    if (verbose){
      cat(c("iter:", j, "; new loglike:", round(loglike.new,4), "; loglike change:", round(loglike.diff,7), "\n"))
    }
    loglike[j]  <- loglike.new
    alpha1.current <- alpha1.new
    alpha2.current <- alpha2.new
    alpha3.current <- alpha3.new
    beta.current   <- beta.new
    loglike.current <- loglike.new
    j <- j+1
  }

  if (!all(loglike[1:(j-1)] == cummax(loglike[1:(j-1)]))){
    cat("Loglikelihood is not strictly monotonically increasing!", "\n")
  }
  noparams <- 4
  AIC <- (-2*loglike.new + noparams * 2)
  BIC <- (-2*loglike.new + noparams * log(n))

  result<-list(estimate = c("alpha1"=alpha1.new,
                            "alpha2"=alpha2.new,
                            "alpha3"=alpha3.new,
                            "beta"  =beta.new),
               loglike    = loglike.new,
               ll         = loglike[1:(j-1)],
               df         = noparams,
               AIC        = AIC,
               BIC        = BIC,
               iter       = (j-1),
               call       = tempcall)
  structure(result, class = c('BGE'))
}


BGE.Estep.function <- function(data, alpha, beta){
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
    #if ( is.nan(s[i]) || is.na(s[i]) ){stop("warings3")}
  }
  return(list("s"=s, "logs"=logs, "logy1s"=logy1s, "logy2s"=logy2s))
}

BGE.data.loglikelihood <- function(data, alpha, beta){
  n <- dim(data)[1]
  den <- rep(0,n)
  for (i in seq_len(n)) {
    den[i] <- dbivgamma(c(data[i,1],data[i,2]) ,
                        alpha=c(alpha[1],alpha[2],alpha[3]),
                        beta=beta,
                        log=TRUE)$value
    if (abs(den[i])==Inf){stop("warings1")}
  }
  return(sum(den))
}

print.BGE <- function (x, ...){ print(x$estimate) }



