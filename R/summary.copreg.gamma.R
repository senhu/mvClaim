#' Summarizing Bivariate Copula Regression with Gamma Margins Model Fits
#'
#' Summary method for class "\code{copreg.gamma}"
#'
#' @param object An object of class "\code{copreg.gamma}" resulting of a call to \code{copreg.gamma}.
#' @param x An object of class "\code{summary.copreg.gamma}", usually a result of a call to \code{summary.copreg.gamma}.
#' @param digits The number of significant digits to use when printing.
#' @param ... Further arguments passed to or from other methods.
#'
#' @examples
#'
#' \donttest{
#' mod <- copreg.gamma(f1 = y1 ~ x1+x2,
#'                     f2 = y2 ~ x1+x2,
#'                     copula =  gumbelCopula(dim=2),
#'                     data = simdat.mcgr)
#' summary(mod)
#' }
#'
#' @export summary.copreg.gamma
#' @export

summary.copreg.gamma <- function(object, ...){
  title <- paste("Bivariate Copula Regression with Gamma Margins")

  obj <- list(title = title,
              copula = object$copula,
              copula.param = object$copula.param,
              n = object$n,
              loglike = object$loglike,
              df = object$df,
              bic = object$BIC,
              aic = object$AIC,
              coefficients = object$coefficients,
              formula = object$formula
  )
  class(obj) <- "summary.copreg.gamma"
  return(obj)
}

#' @rdname summary.copreg.gamma
#' @export print.summary.copreg.gamma
#' @export

print.summary.copreg.gamma <- function(x, digits = getOption("digits"), ...)
{
  txt <- paste(rep("-", min(nchar(x$title), getOption("width"))), collapse = "")
  cat(txt, "\n")
  cat(x$title, "\n")
  cat(txt, "\n")
  #
  cat("\n")
  cat(paste0( class(x$copula), paste0("(",round(x$copula.param,4),") "), "regression:"), "\n")
  cat("\n")
  #
  for (k in 1:2){
    tempcoef <- t(as.matrix(c(x$coefficients$beta[[k]],x$coefficients$shape[k])))
    colnames(tempcoef)[length(tempcoef)]<- "shape"
    rownames(tempcoef) <- names(x$formula)[k]
    print(tempcoef, digits = digits)
  }
  cat("\n")
  #
  tab <- data.frame("log-likelihood" = x$loglike, "n" = x$n,
                    "df" = x$df, "AIC" = x$aic, "BIC" = x$bic,
                    row.names = "", check.names = FALSE)
  print(tab, digits = digits)
  invisible(x)
}
