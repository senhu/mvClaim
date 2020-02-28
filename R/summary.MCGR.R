#' Summarizing Mixture of Copula Regressions with Gamma Marginal Distributions Model Fits
#'
#' Summary method for class "\code{MCGR}"
#'
#' @param object An object of class "\code{MCGR}" resulting of a call to \code{MCGR}.
#' @param x An object of class "\code{summary.MCGR}", usually a result of a call to \code{summary.MCGR}.
#' @param digits The number of significant digits to use when printing.
#' @param ... Further arguments passed to or from other methods.
#'
#' @examples
#' \donttest{
#' mod <- MCGR(copula = list(gumbelCopula(dim=2), frankCopula(dim=2)),
#'             f1 = y1 ~ x1+x2,
#'             f2 = y2 ~ x1+x2,
#'             G = 2,
#'             gating = "C",
#'             data=simdat.mcgr)
#' summary(mod)
#' }
#'
#' @export summary.MCGR
#' @export

summary.MCGR <- function(object, ...){
  title <- paste("Mixture of copula regressions with gamma margins fitted by EM algorithm")
  modelfullname <- paste0("Mixture of ", paste0(sapply(object$copula, class), paste0("(",sapply(object$copula.param,round, 4),")"),collapse = " and "), collapse="")

  gate <- ifelse(object$gating=="V", as.character(object$formula$gating[2]), "None")
  margin1 <- ifelse(is.null(object$formula$f1), "None", as.character(object$formula$f1[3]))
  margin2 <- ifelse(is.null(object$formula$f2), "None", as.character(object$formula$f2[3]))

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
              margin1 = margin1,
              margin2 = margin2,
              classification = object$class
  )
  class(obj) <- "summary.MCGR"
  return(obj)
}

#' @rdname summary.MCGR
#' @export print.summary.MCGR
#' @export

print.summary.MCGR <- function(x, digits = getOption("digits"), ...)
{
  txt <- paste(rep("-", min(nchar(x$title), getOption("width"))), collapse = "")
  cat(txt, "\n")
  cat(x$title, "\n")
  cat(txt, "\n")
  #
  cat("\n")
  cat(paste0(x$fullmodelName," model, ",
             x$G, ifelse(x$G > 1, " components", " component"), ":"), "\n")
  cat("\n")
  #
  cat("Gating Network Covariates:", x$gate, "\n")
  cat("Margin 1 Covariates:", x$margin1, "\n")
  cat("Margin 2 Covariates:", x$margin2, "\n")
  cat("\n")

  tab <- data.frame("log-likelihood" = x$loglike, "n" = x$n,
                    "df" = x$df, "AIC" = x$aic, "BIC" = x$bic,
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
