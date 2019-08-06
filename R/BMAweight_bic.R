#' Approximate model posterior probability using BIC.
#'
#' This function calculates approximately model posterior probability given BIC values, which can then be used for BMA.
#'
#' @param BICvec a vector of all BIC values used to calculate model weight
#' @return the approximated posterior model probability, i.e. model weight used in BMA.
# #' @seealso \code{\link[BMA]{bma.glm}}
#' @examples
#' BMAweight_bic(c(10,11,9,15,9.5))
#'
#' @export BMAweight_bic

BMAweight_bic <- function(BICvec){
  newvec <- BICvec-(min(BICvec))
  expvec <- exp(-0.5 * newvec)
  weight.vec <- expvec / sum(expvec)
  return(weight.vec)
}
