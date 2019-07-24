#' Summarizing Mixture of Bivariate Gamma Clustering Model Fits
#'
#' Summary method for class "\code{MBGC}"
#'
#' @param object An object of class "\code{MBGC}" resulting of a call to \code{MBGC}.
#' @param x An object of class "\code{summary.MBGC}", usually a result of a call to \code{summary.MBGC}.
#' @param digits The number of significant digits to use when printing.
#' @param ... Further arguments passed to or from other methods.
#'
#' @examples
#' \donttest{
#' clust1 <- MBGC(modelName = "CC", y=c("y1","y2"),
#'                G=2, gating = "C", data=gatingsim)
#' summary(clust1)
#' clust2 <- MBGC(modelName = "CI", y=c("y1","y2"),
#'                G=2, gating = "C", data=gatingsim)
#' summary(clust2)
#' clust4 <- MBGC(modelName = "CC", y=c("y1","y2"),
#'                G=2, gating = ~w1+w2+w3, data=gatingsim)
#' summary(clust4)
#' }
#'
#' @export summary.MBGC
#' @export

summary.MBGC <- function(object, ...){
  title <- paste("Mixture of bivariate gamma clustering (MBGC) fitted by EM algorithm")
  modelfullname <- paste0(object$gating, object$modelName, collapse="")

  obj <- list(title = title,
              fullmodelName = modelfullname,
              n = object$n,
              G = object$G,
              loglike = object$loglike,
              df = object$df,
              bic = object$BIC,
              aic = object$AIC,
              pro = object$pro,
              classification = object$class
  )
  class(obj) <- "summary.MBGC"
  return(obj)
}

#' @rdname summary.MBGC
#' @export print.summary.MBGC
#' @export

print.summary.MBGC <- function(x, digits = getOption("digits"), ...)
{
  txt <- paste(rep("-", min(nchar(x$title), getOption("width"))), collapse = "")
  cat(txt, "\n")
  cat(x$title, "\n")
  cat(txt, "\n")
  #
  cat("\n")
  cat(paste0("MBGC ", x$fullmodelName," model with ",
             x$G, ifelse(x$G > 1, " components", " component"), ":"), "\n")
  cat("\n")
  #
  tab <- data.frame("log-likelihood" = x$loglike, "n" = x$n,
                    "df" = x$df, "BIC" = x$bic, "AIC" = x$aic,
                    row.names = "", check.names = FALSE)
  print(tab, digits = digits)
  #
  cat("\nClustering table:")
  print(table(factor(x$classification,
                     levels = { l <- seq_len(x$G)
                     l })),
        digits = digits)
  #
  invisible(x)
}
