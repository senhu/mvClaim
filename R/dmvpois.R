#' Estimate of the density for a given vector lambda of multivariate Poisson distribution (MVP)
#'
#' @param y a vector of MVP sample
#' @param lambda underlying Poisson rates
#' @param method "analytical" or "recursive"
#' @param log whether the calculated value on log scale
#' @param n.mc number of Monte Carlo samples if method='MC"
#'
#' @return Estimate of the probability P(x) for a given vector lambda of Poisson rates
#' @note This function is borrowed from "mvrpois" package available from GitHub.
#' @examples
#' dmvpois(y=c(2,2), lambda=c(1,1,1), method="recursive", log = F)
#' # analytical method only works for d>=3, i.e. not working for bivariate case
#'
#' dmvpois(y=c(2,2,2), lambda=c(1,1,1,1,1,1), method="recursive", log = F)
#' dmvpois(y=c(2,2,2), lambda=c(1,1,1,1,1,1), method="analytical", log = F)
#' dmvpois(y=c(2,2,2), lambda=c(1,1,1,1,1,1), method="MC", log = F)
#'
#' @export


dmvpois <- function(y, lambda, method, log, n.mc = 10000) {
  switch(method,
         MC = {
           A <- mvp.matrix(length(y))
           probability <- sum(colSums((A %*% matrix(rpois(n.mc * length(lambda), lambda),
                                            length(lambda))) == y) == length(y)) / n.mc
           if (log) probability <- log(probability)
         },
         recursive =  {
           if (is.null(dim(y))) {  # x is a single sample
             dict <<- array(-1, dim = y + 1)
             A <- mvp.matrix(length(y))

             if (log) probability <- recur.prob.log(y, A, log(lambda))
             else           probability <- recur.prob(y, A, lambda)
           } else {               # x is an array of samples
             dict <<- array(-1, dim = apply(y, 2, max) + 1 + 1)
             A <- mvp.matrix(ncol(y))

             if (log) probability <-  sum(apply(y, 1, function(x) recur.prob.log(y, A, log(lambda))))
             else           probability <- prod(apply(y, 1, function(x) recur.prob(y, A, lambda)))
           }
         },
         analytical = {
           A <- mvp.matrix(length(y))
           A2 <- A[, (length(y) + 1):ncol(A)]
           tmp <- y * A2 - (A2 - 1) * max(y);

           maximums <- apply(tmp, 2, min)
           y2 <- rep(0, length(maximums))
           probs <- rep(-Inf, prod(maximums + 1))

           j <- 1
           while (!(all(y2 == maximums))) {
             y1 <- y - A2 %*% y2
             if (sum(y1 < 0 ) == 0) probs[j] <- sum(dpois(c(y1, y2), lambda, log = T))

             y2[1] <- y2[1] + 1
             for (i in 1:(length(y2) - 1)) {
               if (y2[i] > maximums[i]) {
                 y2[i] <- 0
                 y2[i + 1] <- y2[i + 1] + 1
               }
             }
             j <- j + 1

           }
           probs[j + 1] <- sum(dpois(c(y - A2 %*% y2, y2), lambda, log = T))
           log.prob <- matrixStats::logSumExp(probs)

           if (log) probability <- log.prob
           else           probability <- exp(log.prob)
         },
         stop('Unknown method. Please use "MC", "recursive", or "analytical"')
  )
  return(probability)
}


############################## PRIVATE HELPER FUNCTIONS ##############################

recur.prob.log <- function(y, A, loglambda) {
  if (sum(y < 0) > 0) return(-Inf)
  if (sum(y) == 0)    return(-exp(matrixStats::logSumExp(loglambda)))

  nnzf <- which.max(y > 0)
  indices <- which(A[nnzf, ] > 0)
  logp <- rep(-Inf, length(indices))
  for (i in 1:length(indices)) {
    idx <- indices[i]
    mvp <- matrix(y - A[, idx] + 1, 1)

    if (length(dict[mvp]) == 0) next
    if (dict[mvp] == -1) {
      result <- recur.prob.log(y - A[, idx], A, loglambda)
      dict[mvp] <<- result
    }
    logp[i] <- loglambda[idx] + dict[mvp]
  }
  return(matrixStats::logSumExp(logp) - log(y[nnzf]))
}

recur.prob <- function(y, A, lambda) {
  if (sum(y < 0) > 0) return(0)
  if (sum(y) == 0) return(exp(-sum(lambda)))

  nnzf <- which.max(y > 0)
  p <- 0
  for (idx in which(A[nnzf, ] > 0)) {
    mvp <- matrix(y - A[, idx] + 1, 1)

    if (length(dict[mvp]) == 0) next
    if (dict[mvp] == -1) {
      result <- recur.prob(y - A[, idx], A, lambda)
      dict[mvp] <<- result
    }
    p <- p + lambda[idx] * dict[mvp]
  }
  return(p / y[nnzf])
}


# mvp.matrix
# m dimension of the m-variate Poisson
# return  the MVP matrix A from X = AY (Karlis 2005), where X is MVP-distributed, and Y is a vector of independent Poissons.
mvp.matrix <- function(m) {
  if (m < 2) {
    stop("Multivariate Poisson is only defined for m > 2. Are you looking for pois, the univariate Poisson distribution?")
  } else {
    combs <- cbind(rbind(1:m, 1:m), combn(1:m, 2))
    k <- ncol(combs)
    A <- matrix(0, m, k)
    for (i in 1:m)
      for (j in 1:k)
        if (combs[1, j] == i | combs[2, j] == i)
          A[i, j] = 1
    return(A)
  }
}

