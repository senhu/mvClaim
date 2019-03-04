mvClaim_icl <- function (object) {
  n <- object$n
  z <- object$z
  if (is.null(z))  z <- matrix(1, nrow = n, ncol = 1)
  C <- matrix(0, n, ncol(z))
  for (i in 1:n) C[i, which.max(z[i, ])] <- 1
  object$BIC + 2 * sum(C * ifelse(z > 0, log(z), 0))
}
