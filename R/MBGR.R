#' Mixture of Bivariate Gamma Clustering/Regressions
#'
#' Mixture of bivariate gamma regressions, or model-based clustering with bivariage gamma distributions and covariates, for various parsimonious model types. Models are estimated by EM algorithm.
#'
#' @param modelName A character string indicating which model to be fitted. Need to be one of "EE", "EI", "IE".
#' @param y A vector of character strings indicating which variables in the data are treated as response or dependent variables.
#' @param data A matrix or data frame of observations. Categorical variables are allowed as covariates.
#' @param G An interger vector specifying the number of mixture components (clusters).
#' @param f1 A regression formula for the \eqn{latex}{\alpha_1} parameter in the bivariate gamma distribution. Note that, depending on the model type, might not be necessary to provide it.
#' @param f2 A regression formula for the \eqn{latex}{\alpha_2} parameter in the bivariate gamma distribution. Note that, depending on the model type, might not be necessary to provide it.
#' @param f3 A regression formula for the \eqn{latex}{\alpha_3} parameter in the bivariate gamma distribution. Note that, depending on the model type, might not be necessary to provide it.
#' @param f4 A regression formula for the \eqn{latex}{\beta} parameter in the bivariate gamma distribution. Note that, depending on the model type, might not be necessary to provide it.
#' @param gating Specifies the gating network in the MoE model, can be "C", "E" or a regression formula.
#' @param initialization Specifies initialization method for EM algorithm. The default is "\code{mclust}".
#' @param maxit A parameter that controls the number of maximum iteration in the EM algorithm. The default is 100.
#' @param tol A parameter that controls the convergence tolerance in the EM algorithm. The default is 1e-5.
#' @param verbose A logical controlling whether estimations in each EM iteration are shown in the fitting process. The default is TRUE.
#'
#' @return An object of class \code{BGR} providing the estimation results.
#'   The details of the output components are:
#'   \item{call}{The matched call.}
#'   \item{expert.coef}{The estimated coefficients in the expert network, if exists.}
#'   \item{gating.coef}{The estimated coefficients in the gating network, if exists.}
#'   \item{alpha1}{The estimated alpha1 values.}
#'   \item{alpha2}{The estimated alpha2 values.}
#'   \item{alpha3}{The estimated alpha3 values.}
#'   \item{beta}{The estimated beta values.}
#'   \item{pro}{A vector whose g-th component is the mixing proportion for the g-th component of the mixture model.}
#'   \item{z}{A matrix whose [i,g]-th entry is the probability that observation i in the data belongs to the g-th group.}
#'   \item{class}{The classification corresponding to z.}
#'   \item{fitted.values}{The fitted values of the regression.}
#'   \item{residuals}{The residuals of the regression}
#'   \item{loglike}{The final estimated maximum log-likelihood value.}
#'   \item{ll}{The sequence of log-likelihood values in the EM algorithm fitting process.}
#'   \item{df}{The number of estimated parameters.}
#'   \item{AIC}{AIC values.}
#'   \item{BIC}{BIC values.}
#'   \item{Hessian}{The Hessian matrix at the estimated values}
#'   \item{iter}{Total iteration numbers.}
#'   \item{formula}{The formulas used in the regression.}
#'   \item{y}{The input response data.}
#'   \item{n}{The number of observations in the data.}
#'   \item{gating.model}{The binomial/multinomial regression model in the gating network.}
#'   \item{Model.Matrix}{The used model matrix for each regression formula.}
#'
#' @examples
#'
#' \dontrun{
#' m1 <- MBGR(modelName = "VV",
#'            y=c("y1","y2"), data = fullsim, G=2,
#'            f1     = ~ w1 + w2,
#'            f2     = ~ w2 + w3,
#'            f3     = ~ w1 + w2 + w3,
#'            f4     = ~ w1 + w2 + w3,
#'            gating = "C")
#' }
#'
#' @importFrom stats glm optim uniroot model.matrix Gamma formula runif coef optimHess
#' @importFrom mclust mclustBIC
#' @export

MBGR <- function(modelName = c("VC","VI","VV","VE", "CV", "IV", "EV", "EC", "CE"),
                 y,
                 data,
                 G,
                 f1,
                 f2,
                 f3,
                 f4,
                 gating,
                 initialization = "mclust",
                 maxit   = 200,
                 tol     = 1e-5,
                 verbose = TRUE){
  switch(modelName,
         VC = {
           if (is.null(f1)){ stop("f1 must be supplied for VC model type") }
           if (is.null(f2)){ stop("f2 must be supplied for VC model type") }
           if (is.null(f3)){ stop("f3 must be supplied for VC model type") }
           res <- MBGR_VC(y=y,data=data,G=G,
                          f1=f1,f2=f2,f3=f3,gating=gating,
                          initialization=initialization,
                          maxit=maxit,tol=tol,verbose=verbose)
         },
         CV = {
           if (is.null(f4)){ stop("f4 must be supplied for CV model type") }
           res <- MBGR_CV(y=y,data=data,G=G,f4=f4,gating=gating,
                          maxit=maxit,tol=tol,
                          initialization=initialization,
                          verbose=verbose)
         },
         VV = {
           if (is.null(f1)){ stop("f1 must be supplied for VV model type") }
           if (is.null(f2)){ stop("f2 must be supplied for VV model type") }
           if (is.null(f3)){ stop("f3 must be supplied for VV model type") }
           if (is.null(f4)){ stop("f4 must be supplied for VV model type") }
           res <- MBGR_VV(y=y,data=data,G=G,
                          f1=f1,f2=f2,f3=f3,f4=f4,gating=gating,
                          initialization=initialization,
                          maxit=maxit,tol=tol,verbose=verbose)
         },
         VI = {
           if (is.null(f1)){ stop("f1 must be supplied for VI model type") }
           if (is.null(f2)){ stop("f2 must be supplied for VI model type") }
           if (is.null(f3)){ stop("f3 must be supplied for VI model type") }
           res <- MBGR_VI(y=y,data=data,G=G,
                          f1=f1,f2=f2,f3=f3,gating=gating,
                          initialization=initialization,
                          maxit=maxit,tol=tol,verbose=verbose)
         },
         IV = {
           if (is.null(f4)){ stop("f4 must be supplied for IV model type") }
           res <- MBGR_IV(y=y,data=data,G=G,f4=f4,gating=gating,
                          maxit=maxit,tol=tol,
                          initialization=initialization,
                          verbose=verbose)
         },
         VE = {
           if (is.null(f1)){ stop("f1 must be supplied for VE model type") }
           if (is.null(f2)){ stop("f2 must be supplied for VE model type") }
           if (is.null(f3)){ stop("f3 must be supplied for VE model type") }
           if (is.null(f4)){ stop("f4 must be supplied for VE model type") }
           res <- MBGR_VE(y=y,data=data,G=G,
                          f1=f1,f2=f2,f3=f3,f4=f4,gating=gating,
                          initialization=initialization,
                          maxit=maxit,tol=tol,verbose=verbose)
         },
         EV = {
           if (is.null(f1)){ stop("f1 must be supplied for EV model type") }
           if (is.null(f2)){ stop("f2 must be supplied for EV model type") }
           if (is.null(f3)){ stop("f3 must be supplied for EV model type") }
           if (is.null(f4)){ stop("f4 must be supplied for EV model type") }
           res <- MBGR_EV(y=y,data=data,G=G,
                          f1=f1,f2=f2,f3=f3,f4=f4,gating=gating,
                          initialization=initialization,
                          maxit=maxit,tol=tol,verbose=verbose)
         },
         EC = {
           if (is.null(f1)){ stop("f1 must be supplied for EC model type") }
           if (is.null(f2)){ stop("f2 must be supplied for EC model type") }
           if (is.null(f3)){ stop("f3 must be supplied for EC model type") }
           res <- MBGR_EC(y=y,data=data,G=G,
                          f1=f1,f2=f2,f3=f3,gating=gating,
                          initialization=initialization,
                          maxit=maxit,tol=tol,verbose=verbose)
         },
         CE = {
           if (is.null(f4)){ stop("f4 must be supplied for CE model type") }
           res <- MBGR_CE(y=y,data=data,G=G,f4=f4,gating=gating,
                          maxit=maxit,tol=tol,
                          initialization=initialization,
                          verbose=verbose)
         },
         stop("invalid model type name") )
  return(res)
}

#' @export
print.MBGR <- function (x, ...){
  modelfullname <- paste0(x$gating, x$modelName, collapse="")
  txt <- paste0("'", class(x)[1], "' model of type '", modelfullname, "' with G = ", x$G)
  cat(txt, "\n")
  cat("\n")
  cat("Available components:\n")
  print(names(x))
  invisible(x)
}

