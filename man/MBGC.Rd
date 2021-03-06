% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MBGC.R
\name{MBGC}
\alias{MBGC}
\alias{MBGC_CC}
\alias{MBGC_CI}
\alias{MBGC_IC}
\title{Mixture of bivariate gamma distributions clustering}
\usage{
MBGC(modelName = c("CC", "CI", "IC"), y, G, gating, data, maxit = 300,
  tol = 1e-05, initialization = "mclust", verbose = FALSE)

MBGC_CC(y, G, gating, data, maxit = 300, tol = 1e-05,
  initialization = "mclust", verbose = FALSE)

MBGC_CI(y, G, gating, data, maxit = 300, tol = 1e-05,
  initialization = "mclust", verbose = FALSE)

MBGC_IC(y, G, gating, data, maxit = 300, tol = 1e-05,
  initialization = "mclust", verbose = FALSE)
}
\arguments{
\item{modelName}{A character string indicating which model to be fitted. Need to be one of "CC", "CI", "IC".}

\item{y}{A vector of character strings indicating which variables in the data are treated as response or dependent variables.}

\item{G}{An integer specifying the numbers of mixture components.}

\item{gating}{Specifies the gating network in the MoE model, can be "C", "E" or a regression formula.}

\item{data}{A matrix or data frame of observations. Categorical variables are allowed as covariates.}

\item{maxit}{A parameter that controls the number of maximum iteration in the EM algorithm. The default is 100.}

\item{tol}{A parameter that controls the convergence tolerance in the EM algorithm. The default is 1e-5.}

\item{initialization}{Specifies initialization method for EM algorithm. The default is "\code{mclust}".}

\item{verbose}{A logical controlling whether estimations in each EM iteration are shown in the fitting process. The default is TRUE.}
}
\value{
An object of class \code{MBGC} providing the estimation results.
  The details of the output components are:
  \item{modelName}{A character string denoting the fitted expert model type.}
  \item{gating}{A character string denoting the fitted gating model type.}
  \item{alpha1}{The estimated alpha1 values.}
  \item{alpha2}{The estimated alpha2 values.}
  \item{alpha3}{The estimated alpha3 values.}
  \item{beta}{The estimated beta values.}
  \item{gating.coef}{The regression coefficients in the gating network, if a regression formula is provided in \code{gating}.}
  \item{pro}{A vector whose g-th component is the mixing proportion for the g-th component of the mixture model.}
  \item{z}{A matrix whose [i,g]-th entry is the probability that observation i in the data belongs to the g-th group.}
  \item{class}{The classification corresponding to z.}
  \item{G}{The number of mixture components.}
  \item{loglike}{The final estimated maximum log-likelihood value.}
  \item{ll}{The sequence of log-likelihood values in the EM algorithm fitting process.}
  \item{df}{Number of estimated parameters.}
  \item{AIC}{AIC values.}
  \item{BIC}{BIC values.}
  \item{iter}{Total iteration numbers.}
  \item{formula}{The formulas used in the regression.}
  \item{gating.model}{The final fitted regression model in gating network.}
  \item{y}{The input response data.}
  \item{n}{The number of observations in the data.}
  \item{Hessian}{The Hessian matrix at the estimated values}
  \item{call}{The matched call.}
}
\description{
Estimation using EM algorithm for mixture of bivariate gamma distributions
}
\examples{
\donttest{
clust1 <- MBGC(modelName = "CC", y=c("y1","y2"),
               G=2, gating = "C", data=gatingsim, verbose=FALSE)
clust1
clust2 <- MBGC(modelName = "CI", y=c("y1","y2"),
               G=2, gating = "C", data=gatingsim, verbose=FALSE)
clust2
clust3 <- MBGC(modelName = "IC", y=c("y1","y2"),
               G=2, gating = "C", data=gatingsim, verbose=FALSE)
clust3
clust4 <- MBGC(modelName = "CC", y=c("y1","y2"),
               G=2, gating = ~w1+w2+w3, data=gatingsim, verbose=FALSE)
clust4
clust5 <- MBGC(modelName = "CI", y=c("y1","y2"),
               G=2, gating = ~w1+w2+w3, data=gatingsim, verbose=FALSE)
clust5
clust6 <- MBGC(modelName = "IC", y=c("y1","y2"),
               G=2, gating = ~w1+w2+w3, data=gatingsim, verbose=FALSE)
clust6
}

}
