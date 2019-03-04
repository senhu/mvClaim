my.gamma.shape <- function(object){
  res <- tryCatch({
    MASS::gamma.shape(object)$alpha
  },
  error=function(e){
    1/(summary(object)$dispersion)
    #print("approximation in gamma shape parameter")
  })
  return(res)
}

# my.gamma.shape <- function (object,
#                             it.lim = 10,
#                             eps.max = .Machine$double.eps^0.25,
#                             verbose = FALSE, ...) {
#   if (is.null(object$y))
#     object <- update(object, y = TRUE)
#   y <- object$y
#   A <- object$prior.weights
#   if (is.null(A))
#     A <- rep(1, length(y))
#   u <- object$fitted.values
#   Dbar <- object$deviance/object$df.residual
#   alpha <- (6 + 2 * Dbar)/(Dbar * (6 + Dbar))
#   if (verbose) {
#     message(gettextf("Initial estimate: %s", format(alpha)),
#             domain = NA)
#     utils::flush.console()
#   }
#   fixed <- -y/u - log(u) + log(A) + 1 + log(y + (y == 0))
#   eps <- 1
#   itr <- 0
#   while (abs(eps) > eps.max && (itr <- itr + 1) <= it.lim) {
#     sc <- sum(A * (fixed + log(alpha) - digamma(A * alpha)))
#     trigamma.temp <- trigamma(A * alpha)
#     if (anyNA(trigamma.temp)){trigamma.temp<- replace(trigamma.temp, is.na(trigamma.temp), 0)}
#     inf <- sum(A * (A * trigamma.temp - 1/alpha))
#     alpha <- alpha + (eps <- sc/inf)
#     if (verbose) {
#       message(gettextf("Iter. %d Alpha: %s", itr, format(alpha)),
#               domain = NA)
#       utils::flush.console()
#     }
#   }
#   if (itr > it.lim)
#     warning("iteration limit reached")
#   res <- list(alpha = alpha, SE = sqrt(1/inf))
#   class(res) <- "gamma.shape"
#   res
# }
