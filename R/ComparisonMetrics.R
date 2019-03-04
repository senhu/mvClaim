require(entropy)
require(epiR)
require(scoringRules)

ComparisonMetrics <- function(actual, pred, exposure=NULL){
  # Gini
  if (is.null(exposure)){
    gini.index <- Gini(actual, pred)
  } else gini.index <- Gini(actual, pred, exposure)

  # K-S test
  #KStest <- ks.test(actual, pred)
  #KStest.statistic <- as.numeric(KStest$statistic)
  #KStest.pvalue <- as.numeric(KStest$p.value)
  # Wasserstein
  Was <- wasserstein1d(actual, pred)
  # concordanace correlation coefficient
  CCC <- as.numeric(epi.ccc(actual, pred)$rho.c[1])

  # KL divergence; empirical KL divergence
  KL <- KL.empirical(actual, pred)
  # freq1 <- freqs(actual, method="ML")
  # freq2 <- freqs(pred, method="ML")
  # KL.2 <- KL.plugin(freq1, freq2)
  # MSE

  MSE <- mean( (actual-pred)^2 )
  RMSE <- sqrt(MSE)

  meanScore <- mean(sapply(actual, crps_sample, dat=pred))

  ResTable <- c("Gini"=gini.index,
                # "KStest.statistic"=KStest.statistic,
                # "KStest.p"=KStest.pvalue,
                "Wasserstein"=Was,
                "CCC"=CCC,
                "KL.divergence"=KL,
                #"MSE"=MSE,
                "RMSE"=RMSE,
                "Scoring"=meanScore)

  return(round(ResTable,4))
}
#AllComparisonMetrics(actual, pred, actual)


wasserstein1d <- function(a, b, p=1, wa=NULL, wb=NULL) {
  m <- length(a)
  n <- length(b)
  stopifnot(m > 0 && n > 0)
  if (m == n && is.null(wa) && is.null(wb)) {
    return(mean(abs(sort(b)-sort(a))^p)^(1/p))
  }
  if (is.null(wa)) {wa <- rep(1,m)}
  if (is.null(wb)) {wb <- rep(1,n)}
  stopifnot(length(wa) == m && length(wb) == n)
  ua <- (wa/sum(wa))[-m]
  ub <- (wb/sum(wb))[-n]
  cua <- c(cumsum(ua))
  cub <- c(cumsum(ub))
  temp <- cut(cub,breaks=c(-Inf,cua,Inf))
  arep <- table(temp) + 1
  temp <- cut(cua,breaks=c(-Inf,cub,Inf))
  brep <- table(temp) + 1
  # we sum over rectangles with cuts on the vertical axis each time one of the two ecdfs makes a jump
  # xrep and yrep tell us how many times each of the x and y data have to be repeated in order to get the points on the horizontal axis
  # note that sum(xrep)+sum(yrep) = m+n-1 (we do not count the height-zero final rectangle where both ecdfs jump to 1)

  aa <- rep(sort(a), times=arep)
  bb <- rep(sort(b), times=brep)

  uu <- sort(c(cua,cub))
  uu0 <- c(0,uu)
  uu1 <- c(uu,1)
  areap <- sum((uu1-uu0)*abs(bb-aa)^p)^(1/p)
  #  print(rbind(uu1-uu0, pmax(aa,bb)-pmin(aa,bb)))
  return(areap)
}


Gini <- function(a, p, exposure = NULL, plot=FALSE) {
  if (length(a) !=  length(p)){stop("Actual and Predicted need to be equal lengths!")}
  if (is.null(exposure)){
    temp.df <- data.frame(actual = a, pred = p, range=c(1:length(a)))

    # re-order the data based on "pred"
    # the only place "pred" used
    temp.df <- temp.df[order(temp.df$pred, temp.df$range, decreasing=TRUE),]

    # accumulated proportion of data (number of how many data points in the data)
    # x-axis
    population.delta <- 1 / length(a)
    null.losses <- rep(population.delta, length(a))
    null.losses.new <- cumsum(null.losses)
    null.losses.new <- c(0, null.losses.new)

    # convert "actual loss" to incremental proportional loss
    # y-axis
    total.losses <- sum(a)
    accum.losses <- temp.df$actual / total.losses
    accum.losses.new <- cumsum(accum.losses)
    accum.losses.new <- c(0, accum.losses.new)

    if (plot){
      plot(null.losses.new, accum.losses.new, type="b", lwd=2, col="blue", main="Lorenz curve",
           xlab="Proportion of data accumulated (sorted by predicitons in decreasing order)",
           ylab="Proportion of losses accumulated")
      abline(a=0,b=1,lwd=2)
    }
    sum.fun          <- rep(0, length(a))
    for (j in c(1:length(a))){
      sum.fun[j]     <- (accum.losses.new[j+1]+accum.losses.new[j])*(null.losses.new[j+1]-null.losses.new[j])/2
    }
    gini             <- sum(sum.fun)
  }
  if (!is.null(exposure)) {
    temp.df          <- data.frame(actual = a, pred = p, exposure = exposure, range=c(1:length(a)))
    # re-order the data based on "pred"
    # the only place "pred" used
    temp.df          <- temp.df[order(temp.df$pred, temp.df$range, decreasing=TRUE),]

    # accumulated proportion of exposure (i.e. data)
    # x-axis
    total.expo       <- sum(exposure)
    prop.expo        <- temp.df$exposure / total.expo
    accum.expo       <- cumsum(prop.expo)
    accum.expo.new   <- c(0, accum.expo)

    # convert "actual loss (or observed)" to incremental proportional loss
    # y-axis
    total.losses     <- sum(a)
    prop.losses      <- temp.df$actual / total.losses
    accum.losses     <- cumsum(prop.losses)
    accum.losses.new <- c(0, accum.losses)

    if (plot==TRUE){
      plot(accum.expo.new, accum.losses.new, type="b", lwd=2, col="blue", main="Lorenz curve",
           xlab="Proportion of exposure accumulated (sorted by predicitons in decreasing order)",
           ylab="Proportion of losses accumulated")
      abline(a=0,b=1,lwd=2)
    }
    sum.fun <- rep(0, length(a))
    for (j in c(1:length(a))){
      sum.fun[j]     <- (accum.losses.new[j+1]+accum.losses.new[j])*(accum.expo.new[j+1]-accum.expo.new[j])/2
    }
    gini             <- sum(sum.fun)
  }
  return(gini)
}
#
# a <- c(1:10)
# p <- c(1,2,2,3,2,5,7,8,9,8)
# exposure <- c(rep(.7,5), rep(1, 5))
# Gini(a, p, exposure, plot=T)
# Gini(a, p, plot=T)


library(transport)

mvComparisonMetric <- function(act, pred, plot=FALSE){
  x <- pp(as.matrix(act))
  y <- pp(as.matrix(pred))
  match <- transport(x,y,p=1)
  if (plot) plot(x,y,match)
  wass <- wasserstein(x,y,p=1,match)
  #wasserstein(x,y,p=1, prob=FALSE)

  return(c(wass))
}


# m4
#
# naviva$COST_AD
#
# m4$fitted.values[,1]
#
# crps_sample(naviva$COST_AD[1], m4$fitted.values[,1])
#
# hist(m4$fitted.values[,1])
# abline(v=naviva$COST_AD[1])
#
#
# mean(sapply(naviva$COST_AD, crps_sample, dat=pred1))
# mean(sapply(naviva$COST_AD, crps_sample, dat=m4$fitted.values[,1]))
# mean(sapply(naviva$COST_AD, crps_sample, dat=m4.2$fitted.values[,1]))
#
#
# mean(sapply(naviva$COST_TPD, crps_sample, dat=pred2))
# mean(sapply(naviva$COST_TPD, crps_sample, dat=m4$fitted.values[,2]))
# mean(sapply(naviva$COST_TPD, crps_sample, dat=m4.2$fitted.values[,2]))
#
#
#
# mean(sapply(naviva$COST_AD, logs_sample, dat=pred1))
# mean(sapply(naviva$COST_AD, logs_sample, dat=m4$fitted.values[,1]))
# mean(sapply(naviva$COST_AD, logs_sample, dat=m4.2$fitted.values[,1]))
#
#
# logs_sample(y = y, dat = sample)
