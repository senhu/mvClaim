#' Mixture of multivariate Poisson regression (mixture of experts)
#'
#' @param K number of mixtures
#' @param l formula for main effect lambda
#' @param lc formula for covariance term lambda
#' @param lp formula for mixing proportion. If NULL, mixing proportion is not regressed against covariates. Default is NULL.
#' @param w weights in GLM (if applicable) such as "exposure". Default is NULL.
#' @param data
#' @param maxit maximum iteration. Default is 300.
#' @param tol tolerance. Default is 1e-5.
#' @param Aitken if TRUE, using Aitken acceleration estimate criterion. Defaul is TRUE.
#' @param verbose if TRUE, print the iteration results. Default is TRUE.
#'
#' @examples
#' cor(cbind(health$doctorco, health$prescrib, health$hospadmi))
#' l <- list(doctorco~sex+age+income,
#'           prescrib~sex+age+income,
#'           hospadmi~sex+age+income)
#' lc <- list(x12~sex+age+income,
#'            x13~sex+age+income,
#'            x23~sex+age+income)
#' lp <- ~sex+age+income
#' data <- health

#---------------------
# to do list
# fix fitted value part
# change this only to 3-dim
#----------------------

MMPR<- function(K=2,
                l,
                lc=NULL,
                lp=NULL,
                w=NULL,
                data,
                maxit=300,
                tol=1e-5,
                Aitken = TRUE,
                verbose=TRUE){

  options(warn=-1) # sets the handling of warning messages. If warn is negative all warnings are ignored. Default is 0

  # definition of function call
  templist<-list( l=l,
                  lc=lc,
                  lp=lp,
                  w = w,
                  data=substitute(data),
                  maxit=maxit,
                  tol=tol,
                  Aitken=Aitken,
                  verbose=verbose)
  tempcall<-as.call( c(expression(MBPR), templist))
  rm(templist)

  m <- length(l); m.c <- length(lc)
  namey <- rep(NULL, m)
  datay <- NULL
  for (i in seq_len(m)){
    namey[i] <- as.character(l[[i]][2])
    datay <- cbind(datay, data[,names(data)==namey[i]])
    colnames(datay)[i] <- paste("y", i, sep="")
  }
  n<-dim(datay)[1]
  p<-dim(as.data.frame(data))[2]
  if (!is.null(w)) {weight.var <- data[, names(data)==w]} else {weight.var <- rep(1, n)}
  data$weight.var <- weight.var
  datam<-data # internal data for estimating main lambda in M-step
  datam<-datam[ , -match(namey,names(datam))] # removing y
  datac<-datam # internal data for estimating covaraince lambda in M-step
  datap<-datam # internal data for estimating mixing proportion in M-step

  # prepare model formula
  l.new <- l
  for (i in seq_len(m)){
    if (as.character(l[[i]][3])=="."){
      l.new[[i]]<-formula(paste(paste("x",i,sep=""),paste(names(data1),'',collapse='+',sep=''),sep='~'))}
    else {l.new[[i]] <- formula(paste(paste("x",i,sep=""),as.character(lc[[i]])[3],sep="~"))}
  }
  if (!is.null(lc)){
    lc.new <- lc
    for (i in seq_len(length(lc))){
      if (as.character(lc[[i]][3])=="."){
        lc.new[[i]]<-formula(paste(as.character(lc[[i]][2]) ,paste(names(datac),"",collapse="+",sep=""),sep="~"))}
      else {lc.new[[i]]<-formula(paste(as.character(lc[[i]][2]),as.character(lc[[i]])[3],sep="~"))}
    }
  }
  if (!is.null(lp)){
    if (as.character(lp[2])=="."){lp.new<-formula(paste( "PCY", paste( names(data1), "", collapse = "+", sep = ""), sep="~"))}
    else {lp.new <- formula(paste("PCY", as.character(lp)[2], sep="~"))}
  }

  # initial values
  for (i in seq_len(length(lc))){
    assign(as.character(lc[[i]][2]),  matrix(0,nrow=n, ncol=K))
  }
  #like<-matrix(0, nrow=n, ncol=K) # likelihood, vector of length n
  #zero <- (rowSums(datay)==0) # obs who has 0 in either dimention
  zero12 <- (datay[,1] == 0 | datay[,2]==0)
  zero13 <- (datay[,1] == 0 | datay[,3]==0)
  zero23 <- (datay[,2] == 0 | datay[,3]==0)
  #p.c <- rep(1/K, K) # initial value for P(C_k)
  tm <- runif(K); p.c <- tm/sum(tm)
  lambda.m <- NULL
  for (i in seq_len(m)){
    lambda.temp <- glm(l[[i]], weights = weight.var, family=poisson, data=data)$fitted
    #lambda.m <- append(lambda.m, list(matrix(rep(lambda.temp, each=K), ncol=K, byrow=TRUE)))
    lambda.m <- append(lambda.m, list(lambda.temp * t(rmultinom(n, size=100, prob = rep(1/K, K))/100)))
  }
  if (is.null(lc)) {
    lambda.c <- rep(list(rep(0,n)),3)
  } else {
    lambda.c <- NULL
    for (i in seq_len(m.c)){
      lambda.c.temp <- rep( max(0.1, min(cor(datay))), n)
      #lambda.c <- append(lambda.c, list(matrix(rep(lambda.c.temp, each=K), ncol=K, byrow=TRUE)))
      lambda.c <- append(lambda.c, list(lambda.c.temp * t(rmultinom(n, size=100, prob = rep(1/K, K))/100)))
    }
  }
  loglike.diff <- 1000.0 # a starting value of likelihood difference
  loglike.current <- 1000.0 # current log likelihood value at current iteration
  i<-0
  # fitting the Double Poisson Model
  if (is.null(lc)) {
    stop("lc is NULL")
  }
  # fitting MP model
  if (!is.null(lc)){

    loglike<-rep(0,maxit)
    loglike.store <- rep(-Inf, 3) # for Aitken
    betap <- NULL

    lambda.m.track.11 <- NULL
    lambda.m.track.12 <- NULL
    lambda.m.track.21 <- NULL
    lambda.m.track.22 <- NULL
    lambda.m.track.31 <- NULL
    lambda.m.track.32 <- NULL
    lambda.c.track.11 <- NULL
    lambda.c.track.12 <- NULL
    lambda.c.track.21 <- NULL
    lambda.c.track.22 <- NULL
    lambda.c.track.31 <- NULL
    lambda.c.track.32 <- NULL
    betam.track <- NULL
    betac.track <- NULL

    while ( (loglike.diff > tol) && (i <= maxit) ) {
      i<-i+1
      loglikelihood <- rep(0, n)
      PCY <- matrix(0, nrow=n, ncol=K)
      # E step
      for (j in 1:n) { # each observation

        #if (!zero[j]) {lbp1 <- rep(0, K); lbp2 <- rep(0, K)}
        #if (zero[j]) {lbp2 <- rep(0, K)}
        ll <- rep(0, K)
        for (k in 1:K) { # for each cluster

          ll[k] <- dmvpois(c(datay[,1][j], datay[,2][j], datay[,3][j]),
                           lambda=c(lambda.m[[1]][j, k],
                                    lambda.m[[2]][j, k],
                                    lambda.m[[3]][j, k],
                                    lambda.c[[1]][j, k],
                                    lambda.c[[2]][j, k],
                                    lambda.c[[3]][j, k]),
                           method="recursive",
                           log=TRUE)

          if (zero12[j]) { x12[j,k] <- 0.0 } else {
            lmp12 <- dmvpois(c(datay[,1][j]-1, datay[,2][j]-1, datay[,3][j]),
                             lambda=c(lambda.m[[1]][j, k],
                                       lambda.m[[2]][j, k],
                                       lambda.m[[3]][j, k],
                                       lambda.c[[1]][j, k],
                                       lambda.c[[2]][j, k],
                                       lambda.c[[3]][j, k]),
                             method="recursive",
                             log=TRUE)
            x12[j,k]<-lambda.c[[1]][j, k]*exp(lmp12) / exp(ll[k])
            #if (is.nan(s[j, k]) && j>=2 && !is.nan(s[j-1, k])){s[j, k]=s[j-1, k]}
            #else if (is.nan(s[j, k])) {s[j, k] = 0}
            if (is.nan(x12[j,k])){x12[j,k]=0}
          }

          if (zero13[j]) { x13[j,k] <- 0.0 } else {
            lmp13 <- dmvpois(c(datay[,1][j]-1, datay[,2][j], datay[,3][j]-1),
                             lambda=c(lambda.m[[1]][j, k],
                                      lambda.m[[2]][j, k],
                                      lambda.m[[3]][j, k],
                                      lambda.c[[1]][j, k],
                                      lambda.c[[2]][j, k],
                                      lambda.c[[3]][j, k]),
                             method="recursive",
                             log=TRUE)
            x13[j,k]<-lambda.c[[2]][j, k]*exp(lmp13) / exp(ll[k])
            #if (is.nan(s[j, k]) && j>=2 && !is.nan(s[j-1, k])){s[j, k]=s[j-1, k]}
            #else if (is.nan(s[j, k])) {s[j, k] = 0}
            if (is.nan(x13[j,k])){x13[j,k]=0}
          }

          if (zero23[j]) { x23[j,k] <- 0.0 } else {
            lmp23 <- dmvpois(c(datay[,1][j], datay[,2][j]-1, datay[,3][j]-1),
                             lambda=c(lambda.m[[1]][j, k],
                                      lambda.m[[2]][j, k],
                                      lambda.m[[3]][j, k],
                                      lambda.c[[1]][j, k],
                                      lambda.c[[2]][j, k],
                                      lambda.c[[3]][j, k]),
                             method="recursive",
                             log=TRUE)
            x23[j,k]<-lambda.c[[3]][j, k]*exp(lmp23) / exp(ll[k])
            #if (is.nan(s[j, k]) && j>=2 && !is.nan(s[j-1, k])){s[j, k]=s[j-1, k]}
            #else if (is.nan(s[j, k])) {s[j, k] = 0}
            if (is.nan(x23[j,k])){x23[j,k]=0}
          }

        } # end of k
        loglikelihood[j] <- log(p.c %*% exp(ll))
        PCY[j,] <- (p.c * exp(ll)) /  (p.c %*% exp(ll))
      } # end of j
      ##### end of E step  ######

      #####   M step  ######
      if (is.null(lp)){ p.c <- colSums(PCY)/n } else {
        mp <- multinom(formula = lp.new, data = datap)
        betap <- append(betap, list(coef(mp)))
        p.c <- colSums(mp$fitted.values)/n
      }

      betam <- rep(list(NULL),m)
      betac <- rep(list(NULL),m)
      for (k in c(1:K)){
        x1<-abs(datay[,1]-x12[,k]-x13[,k])
        x2<-abs(datay[,2]-x12[,k]-x23[,k])
        x3<-abs(datay[,3]-x13[,k]-x23[,k])
        datam$x1 <- x1
        datam$x2 <- x2
        datam$x3 <- x3
        for (h in c(1:m)){
          mod.m<-glm(l.new[[h]], family=poisson, data=datam, weights = PCY[,k], offset = weight.var)
          betam[[h]] <- cbind(betam[[h]], mod.m$coef)
          lambda.m[[h]][,k]<-mod.m$fitted
        }

        # fit model on lambda3
        datac$x12 <- x12[, k]
        datac$x13 <- x13[, k]
        datac$x23 <- x23[, k]
        for (h in c(1:m.c)){
          mod.c<-glm(lc.new[[h]], family=poisson, data=datac, weights = PCY[, k], offset = weight.var)
          betac[[h]] <- cbind(betac[[h]], mod.c$coef)
          lambda.c[[h]][,k]<-mod.c$fitted
        }
      }
      #####   end of M step  ######

      loglike.new <- sum(loglikelihood)
      loglike[i] <- loglike.new
      loglike.store <- c(loglike.store[-1], loglike.new)
      if (Aitken){
        loglike.diff <- abs( upclass::Aitken(loglike.store)$linf - loglike.current )
      } else {
        loglike.diff <- abs( loglike.new - loglike.current )
      }
      loglike.current<-loglike.new

      #	detailed or compressed printing during the EM iterations
      if (verbose) {
        printvector<-c(i, loglike.new, loglike.diff )
        names(printvector)<-c( 'iter', 'loglike', 'Rel.Dif.loglike')}
      print.default( printvector, digits=4 )

      lambda.m.track.11 <- cbind(lambda.m.track.11, lambda.m[[1]][,1])
      lambda.m.track.12 <- cbind(lambda.m.track.12, lambda.m[[1]][,2])
      lambda.m.track.21 <- cbind(lambda.m.track.21, lambda.m[[2]][,1])
      lambda.m.track.22 <- cbind(lambda.m.track.22, lambda.m[[2]][,2])
      lambda.m.track.31 <- cbind(lambda.m.track.31, lambda.m[[3]][,1])
      lambda.m.track.32 <- cbind(lambda.m.track.32, lambda.m[[3]][,2])

      lambda.c.track.11 <- cbind(lambda.c.track.11, lambda.c[[1]][,1])
      lambda.c.track.12 <- cbind(lambda.c.track.12, lambda.c[[1]][,2])
      lambda.c.track.21 <- cbind(lambda.c.track.21, lambda.c[[2]][,1])
      lambda.c.track.22 <- cbind(lambda.c.track.22, lambda.c[[2]][,2])
      lambda.c.track.31 <- cbind(lambda.c.track.31, lambda.c[[3]][,1])
      lambda.c.track.32 <- cbind(lambda.c.track.32, lambda.c[[3]][,2])

      betam.track <- cbind(betam.track, betam[[1]][,1])
      betac.track <- cbind(betac.track, betac[[1]][,1])
    } # loop over EM steps

    #	calculation of BIC and AIC for multivariate Poisson model
    if (is.null(lp)){
      noparams<- K*(mod.m$rank * m + mod.c$rank * m.c) + (K-1)
    } else {
      noparams<- K*(mod.m$rank * m + mod.c$rank * m.c) + mp$rank
    }
    AIC<- -2*loglike[i] + noparams * 2
    BIC<- -2*loglike[i] + noparams * log(n)

    # fitted values
    y1.fitted <- t(p.c %*% t(lambda.m[[1]] + lambda.c[[1]] + lambda.c[[2]]))
    y2.fitted <- t(p.c %*% t(lambda.m[[2]] + lambda.c[[1]] + lambda.c[[3]]))
    y3.fitted <- t(p.c %*% t(lambda.m[[3]] + lambda.c[[2]] + lambda.c[[3]]))

    #  Calculation of output
    result<-list(coefficients=c(list(betam), list(betac), list(betap)),
                 betam=betam,
                 betac=betac,
                 betap=betap,
                 lambda.m=lambda.m,
                 lambda.c=lambda.c,
                 MixProb = p.c,
                 PCY = PCY,
                 fitted.values=data.frame(y1=y1.fitted,
                                          y2=y2.fitted,
                                          y3=y3.fitted),
                 residuals=data.frame(y1=(datay[,1]-y1.fitted),
                                      y2=(datay[,2]-y2.fitted),
                                      y3=(datay[,3]-y3.fitted)),
                 loglikelihood=loglike[i],
                 llseq=loglike[1:i],
                 parameters.number=noparams,
                 AIC=AIC,
                 BIC=BIC,
                 iterations=i,
                 call=tempcall )

  } # end of multivariate BivPois regression model
  options(warn=0)
  class(result)<-c('MMPR', 'glm')
  return(result)
} # end of function
