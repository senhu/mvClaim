#' Summarizing Bivariate Gamma Distribution Fits
#'
#' Summary method for class "\code{BGE}"
#'
#' @param object An object of class "\code{BGE}" resulting of a call to \code{BGE}.
#' @param x An object of class "\code{summary.BGE}", usually a result of a call to \code{summary.BGE}.
#' @param digits The number of significant digits to use when printing.
#' @param ... Further arguments passed to or from other methods.
#'
#' @examples
#' \donttest{
#' dat <- rbivgamma(1000, alpha = c(1,2,0.5), beta=0.1)
#' mod <- BGE(data = dat, verbose = FALSE)
#' summary(mod)
#' }
#'
#' @export summary.BGE
#' @export

summary.BGE <- function(object, ...){
  title <- paste("Bivariate gamma distribution estimated by EM algorithm")
  modelfullname <- "II"
  comp <- object$estimate
  names(comp) <- c("alpha1", "alpha2", "alpha3", "beta")

  obj <- list(title = title,
              fullmodelName = modelfullname,
              n = object$n,
              loglike = object$loglike,
              df = object$df,
              bic = object$BIC,
              aic = object$AIC,
              parameters = comp
  )
  class(obj) <- "summary.BGE"
  return(obj)
}

#' @rdname summary.BGE
#' @export print.summary.BGE
#' @export

print.summary.BGE <- function(x, digits = getOption("digits"), ...)
{
  txt <- paste(rep("-", min(nchar(x$title), getOption("width"))), collapse = "")
  cat(txt, "\n")
  cat(x$title, "\n")
  cat(txt, "\n")
  #
  cat("\n")
  cat("Estimated parameters:", "\n")
  print(x$parameters, digits = digits)
  cat("\n")
  #
  tab <- data.frame("log-likelihood" = x$loglike, "n" = x$n,
                    "df" = x$df, "AIC" = x$aic, "BIC" = x$bic,
                    row.names = "", check.names = FALSE)
  print(tab, digits = digits)
  invisible(x)
}
