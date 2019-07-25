#' Summarizing Mixture of Bivariate Gamma Regressions Model Fits
#'
#' Summary method for class "\code{MBGR}"
#'
#' @param object An object of class "\code{MBGR}" resulting of a call to \code{MBGR}.
#' @param x An object of class "\code{summary.MBGR}", usually a result of a call to \code{summary.MBGR}.
#' @param digits The number of significant digits to use when printing.
#' @param ... Further arguments passed to or from other methods.
#'
#' @examples
#' \donttest{
#' mod <- MBGR(modelName = "CE",
#'             y=c("y1","y2"),
#'             data   = fullsim,
#'             G=2,
#'             f1     = ~ w1 + w2,
#'             f2     = ~ w2 + w3,
#'             f3     = ~ w1 + w2 + w3,
#'             f4     = ~ w1 + w2 + w3,
#'             gating = "C")
#' summary(mod)
#' }
#'
#' @export summary.MBGR
#' @export

summary.MBGR <- function(object, ...){
  title <- paste("Mixture of bivariate gamma regressions fitted by EM algorithm")
  modelfullname <- paste0(object$gating, object$modelName, collapse="")

  gate <- ifelse(object$gating=="V", as.character(object$formula$gating[2]), "None")
  expert1 <- ifelse(is.null(object$formula$f1), "None", as.character(object$formula$f1[2]))
  expert2 <- ifelse(is.null(object$formula$f2), "None", as.character(object$formula$f2[2]))
  expert3 <- ifelse(is.null(object$formula$f3), "None", as.character(object$formula$f3[2]))
  expert4 <- ifelse(is.null(object$formula$f4), "None", as.character(object$formula$f4[2]))

  obj <- list(title = title,
              fullmodelName = modelfullname,
              n = object$n,
              G = object$G,
              loglike = object$loglike,
              df = object$df,
              bic = object$BIC,
              aic = object$AIC,
              pro = object$pro,
              gate= gate,
              expert1 = expert1,
              expert2 = expert2,
              expert3 = expert3,
              expert4 = expert4,
              classification = object$class
  )
  class(obj) <- "summary.MBGR"
  return(obj)
}

#' @rdname summary.MBGR
#' @export print.summary.MBGR
#' @export

print.summary.MBGR <- function(x, digits = getOption("digits"), ...)
{
  txt <- paste(rep("-", min(nchar(x$title), getOption("width"))), collapse = "")
  cat(txt, "\n")
  cat(x$title, "\n")
  cat(txt, "\n")
  #
  cat("\n")
  cat(paste0("MBGR ", x$fullmodelName," model with ",
             x$G, ifelse(x$G > 1, " components", " component"), ":"), "\n")
  cat("\n")
  #
  cat("Gating Network Covariates:", x$gate, "\n")
  cat("Alpha1 Expert Network Covariates:", x$expert1, "\n")
  cat("Alpha2 Expert Network Covariates:", x$expert2, "\n")
  cat("Alpha3 Expert Network Covariates:", x$expert3, "\n")
  cat("Beta Expert Network Covariates:", x$expert4, "\n")
  cat("\n")

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
