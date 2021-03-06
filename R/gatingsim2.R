#' Simulated artificial data with covariates in the gating network.
#'
#' An artificial data set containing the bivariate response \code{y1}, \code{y2}, three covariates \code{w1}, \code{w2}, \code{w3}, and the true label.
#'
#' @format A data frame with 500 rows and 6 variables:
#' \describe{
#'   \item{\code{y1},\code{y2}}{The bivariate response variable.}
#'   \item{\code{w1},\code{w2},\code{w3}}{Three covariates.}
#'   \item{\code{label}}{Their true labels.}
#' }
#' @details This data set of size 500 is artificially generated as following: first \code{w1}, \code{w2}, and \code{w3} are simulated from \code{Normal(0, 1)}.
#' The true number of component is \eqn{G = 2}.
#' The gating is simulated from a logistic regression where the regression coeffcients are
#' \deqn{ logit(p) = 0.4 + 0.7 w1 + 0.5 w2 + 0.6 w3 ,}
#' based on which a binomial simulation is used.
#' The two components were simulated from two bivariate gamma distributions:
#' for component 1:
#' \deqn{ y ~ BG(\alpha1 = 1, \alpha2 = 8, \alpha3 = 5, \beta = 2.5) ,}
#' for component 2:
#' \deqn{ y ~ BG(\alpha1 = 3, \alpha2 = 3, \alpha3 = 1.5, \beta = 1) . }
#' The simulated data set consists of 222 observations from component 1 and 278 observations from component 2.
#'
#' @source Hu, S., Murphy, T. B. and O'Hagan, A. (2019) Bivariate Gamma Mixture of Experts Models for Joint Claims Modeling. To appear.
#'
"gatingsim2"
