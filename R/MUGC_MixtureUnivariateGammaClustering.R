#' Mixture of univariate Gamma clustering
#'
#' @param G number of mixtures
#' @param data
#' @param maxit maximum iteration. Default is 300.
#' @param tol tolerance. Default is 1e-5.
#' @param Aitken if TRUE, using Aitken acceleration estimate criterion. Defaul is TRUE.
#' @param verbose if TRUE, print the iteration results. Default is TRUE.
#'
#'
#'
#------------------------------
# mixtools::gammamixEM # similar in R package
#-----------------------------

UnivariateGammaMixture <- function(data,
                 G,
                 maxit=300,
                 tol=1e-5,
                 Mstep.method,
                 Aitken=FALSE,
                 verbose=TRUE){
  
  options(warn=0)
  n<-length(data)
  p.z.temp <- rep(1/(G+1), (G-1))
  p.z <-  c(p.z.temp, 1-sum(p.z.temp))
  alpha.current <- rep(1, G)
  beta.current <- rep(1,G)
  
  j<-1
  loglike.diff    <- 1000
  loglike.current <- -sqrt(.Machine$double.xmax)
  loglike         <- rep(0,maxit)
  loglike.store   <- rep(-sqrt(.Machine$double.xmax), 3)
  alpha.track <- alpha.current
  beta.track   <- beta.current
  
  while ( (loglike.diff >= tol) && (j <= maxit) ) {
    #-------
    #E step
    #-------
    den <- matrix(0, nrow=n, ncol=G)
    for (i in 1:n){
      for (g in 1:G){
        den[i,g] <- dgamma(data[i], shape=alpha.current[g], rate=beta.current[g]) * p.z[g]
      }
    }
    likelihood <- rowSums(den)
    pzy <- matrix(0, nrow=n, ncol=G)
    for (w in 1:n){  pzy[w,] <- den[w,]/likelihood[w]  }
    
    loglike.new <- sum(log(likelihood))
    loglike[j]      <- loglike.new
    loglike.store <- c(loglike.store[-1], loglike.new)
    if (Aitken){
      loglike.diff <- abs(MoEClust::MoE_aitken(loglike.store)$linf - loglike.current)
    } else {
      loglike.diff <- abs(loglike.new - loglike.current) / ( 1+abs(loglike.new) )
    }
    loglike.current <- loglike.new ; rm(loglike.new)
    
    #-------
    #M step
    #-------
    alpha.new <- rep(0, G)
    beta.new   <- rep(0, G)
    p.z        <- colSums(pzy)/n
    if (Mstep.method=="numerical.uniroot"){
      for (g in seq_len(G)){
        alpha.rootfun <- function(alpha.var, g){
          log(beta.current[g])*sum(pzy[,g]) - digamma(alpha.var)*sum(pzy[,g]) + sum(pzy[,g]*log(data))
        }
        alpha.new[g] <- uniroot(alpha.rootfun, g=g,
                            lower = .Machine$double.eps , 
                            upper = 100)$root
        beta.new[g] <- alpha.new[g]*sum(pzy[,g]) / sum(pzy[,g]*data)
      }
    }
    if (Mstep.method=="numerical.optim.all"){
      for (g in seq_len(G)){
        target.fun <- function(param, who){
          alpha.var <- param[1]
          beta.var <- param[2]
          log(beta.var)*alpha.var*sum(pzy[,who]) - lgamma(alpha.var)*sum(pzy[,who]) + (alpha.var-1)*(pzy[,who]%*%log(data)) - beta.var*(pzy[,who]%*%data)
        }
        target.gr <- function(param, who){
          alpha.var <- param[1]
          beta.var <- param[2]
          c(log(beta.var)*sum(pzy[,who]) - digamma(alpha.var)*sum(pzy[,who]) + (pzy[,who]%*%log(data)),
            alpha.var/beta.var * sum(pzy[,who]) - pzy[,g]%*%data) 
        }
        optim.res <- optim(c(alpha.current[g], beta.current[g]), target.fun, target.gr, 
                           who=g, method = "L-BFGS-B",
                           lower=rep(sqrt(.Machine$double.xmin),2),
                           control=list(fnscale=-1,ndeps=1e-10,maxit=1000000,factr=.Machine$double.xmin))
        alpha.new[g] <- optim.res$par[1]
        beta.new[g] <- optim.res$par[2]
      }
    }
    if (verbose){
      cat(c("iter:", j, "; current loglike:", round(loglike.current,6), "; loglike change:", round(loglike.diff,6), "\n"))
    }
    
    alpha.track   <- rbind(alpha.track, alpha.new)
    beta.track     <- rbind(beta.track, beta.new)
    
    alpha.current <- alpha.new; rm(alpha.new)
    beta.current   <- beta.new; rm(beta.new)
    j<-j+1
  }
  
  noparams <- 5*G - 1
  AIC<- -2*loglike[j-1] + noparams * 2
  BIC<- -2*loglike[j-1] + noparams * log(n)
  
  classification <- apply(pzy, 1, which.max)
  
  #  Calculation of output
  result<-list(alpha=alpha.current,
               beta=beta.current,
               clusterProb = p.z,
               pzy = pzy,
               loglikelihood=loglike[j-1],
               loglike.seq=loglike[1:(j-1)],
               parameters.number=noparams,
               AIC=AIC,
               BIC=BIC,
               classification=classification,
               alpha.trace=alpha.track,
               beta.trace=beta.track,
               iterations=(j-1))
  options(warn=0)
  class(result)<-c('MBGC')
  return(result)
}

