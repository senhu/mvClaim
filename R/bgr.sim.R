#' Simulated artificial data based on bivariate gamma regression (no mixtures).
#'
#' An artificial data set of size 1000 containing the bivariate response \code{y1}, \code{y2}, three covariates \code{w1}, \code{w2}, \code{w3}.
#'
#' @format A data frame with 1000 rows and 5 variables:
#' \describe{
#'   \item{\code{y1},\code{y2}}{The bivariate response variable.}
#'   \item{\code{w1},\code{w2},\code{w3}}{Three covariates.}
#' }
#' @details This data set of size 1000 is artificially generated as following: first \code{w1}, \code{w2}, and \code{w3} are simulated from \code{Normal(0, 1)}.
#' Then the \eqn{\alpha1, \alpha2, \alpha3, \beta} are generated via
#' \deqn{ \alpha1 = exp(1 + 0.1 w1 + 0.3 w2) }
#' \deqn{ \alpha2 = exp(0.1 + 0.1 w2 + 0.1 w3) }
#' \deqn{ \alpha3 = exp(0.5 + 0.1 w1 + 0.2 w2 + 0.4 w3) }
#' \deqn{ \beta = exp(0.2 + 0.1 w1 + 0.1 w2 + 0.2 w3) }
#' For each set of parameters \eqn{(\alpha1, \alpha2, \alpha3, \beta)} a random sample from this bivariate gamma distribution is simulated, i.e. \eqn{y ~ BG(\alpha1, \alpha2, \alpha3, \beta)}.
#'
#' @source Hu, S., Murphy, T. B. and O'Hagan, A. (2019) Bivariate Gamma Mixture of Experts Models for Joint Claims Modeling. To appear.
"bgr.sim"
