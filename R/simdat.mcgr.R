#' Simulated artificial data with covariates for `MCGR`.
#'
#' An artificial data set of size 1000 containing the bivariate response \code{y1}, \code{y2}, two covariates \code{x1}, \code{x2}, and the true label.
#'
#' @format A data frame with 1000 rows and 5 variables:
#' \describe{
#'   \item{\code{y1},\code{y2}}{The bivariate response variable.}
#'   \item{\code{x1},\code{x2}}{Two covariates.}
#'   \item{\code{label}}{Their true labels.}
#' }
#' @details This data set of size 1000 is artificially generated as following: first \code{x1} and \code{x2} are simulated from \code{Normal(0, 1)}.
#' The true number of component is \eqn{G = 2},
#' and the two components are equally mixed, i.e. \eqn{\tau = (0.5, 0.5)}.
#' Component 1 follows a Gumbel(2) copula, while component 2 follows a Frank(3) copula.
#' For the marginal distributions, they all have gamma distributions, with shape parameters being \eqn{\gamma1 = 6, \gamma2 = 3} for component 1, and \eqn{\gamma1 = 2, \gamma2 = 2} for component 2.
#' For component 1:
#' \deqn{ y1 = 1 + 0.1 x1 + 0.1 x2 , }
#' \deqn{ y2 = 2 + 0.2 x1 + 0.2 x2 . }
#' For component 2:
#' \deqn{ y1 = 0.8 + 0.2 x1 + 0.2 x2 , }
#' \deqn{ y2 = 0.8 + 0.2 x1 + 0.2 x2 . }
#' The final simulated data consists of 500 samples from component 1, 500 samples from component 2
#'
#' @source Hu, S. and O'Hagan, A. (2019) Copula Averaging for Tail Dependence in General Insurance Claims Data. To appear.
"simdat.mcgr"
