#' Summarizing Bivariate Gamma Regression Model Fits
#'
#' Summary method for class "\code{BGR}"
#'
#' @param object An object of class "\code{BGR}" resulting of a call to \code{BGR}.
#' @param x An object of class "\code{summary.BGR}", usually a result of a call to \code{summary.BGR}.
#' @param digits The number of significant digits to use when printing.
#' @param ... Further arguments passed to or from other methods.
#'
#' @examples
#'
#' \donttest{
#' mod <- BGR(modelName = "EE",
#'            y=c("y1","y2"),
#'            data   = fullsim,
#'            f1     = ~ w1 + w2,
#'            f2     = ~ w2 + w3,
#'            f3     = ~ w1 + w2 + w3,
#'            f4     = ~ w1 + w2 + w3,
#'            verbose= FALSE)
#' summary(mod)
#' }
#'
#' @export


summary.BGR <- function(object, ...){
  title <- paste("Bivariate gamma regression (BGR) fitted by EM algorithm")

  obj <- list(title = title,
              fullmodelName = object$modelName,
              n = object$n,
              loglike = object$loglike,
              df = object$df,
              bic = object$BIC,
              aic = object$AIC
  )
  class(obj) <- "summary.BGR"
  return(obj)
}

#' @rdname summary.BGR
#' @export

print.summary.BGR <- function(x, digits = getOption("digits"), ...)
{
  txt <- paste(rep("-", min(nchar(x$title), getOption("width"))), collapse = "")
  cat(txt, "\n")
  cat(x$title, "\n")
  cat(txt, "\n")
  #
  cat("\n")
  cat(paste0("BGR ", x$fullmodelName," model :"), "\n")
  cat("\n")
  #
  tab <- data.frame("log-likelihood" = x$loglike, "n" = x$n,
                    "df" = x$df, "BIC" = x$bic, "AIC" = x$aic,
                    row.names = "", check.names = FALSE)
  print(tab, digits = digits)
  invisible(x)
}
