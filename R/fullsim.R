#' Simulated artificial data with covariates in both the gating and expert networks.
#'
#' An artificial data set of size 500 containing the bivariate response \code{y1}, \code{y2}, three covariates \code{w1}, \code{w2}, \code{w3}, and the true labels.
#'
#' @format A data frame with 500 rows and 6 variables:
#' \describe{
#'   \item{\code{y1},\code{y2}}{The bivariate response variable.}
#'   \item{\code{w1},\code{w2},\code{w3}}{Three covariates.}
#'   \item{\code{label}}{Their true labels.}
#' }
#' @details This data set of size 500 is artificially generated as following: first \code{w1}, \code{w2}, and \code{w3} are simulated from \code{Normal(0, 0.3)}.
#' The true number of component is \eqn{G = 2}.
#' For the first component, the parameters are generated as:
#' \deqn{ \alpha1 = exp(1 + 0.2 w1 + 0.2 w2) }
#' \deqn{ \alpha2 = exp(0.1 + 0.1 w2 + 0.1 w3) }
#' \deqn{ \alpha3 = exp(0.5 + 0.2 w1 + 0.2 w2 + 0.2 w3) }
#' \deqn{ \beta = exp(0.2 + 0.1 w1 + 0.1 w2 + 0.2 w3) }
#' For the second component, the parameters are generated as:
#' \deqn{ \alpha1 = exp(0.1 + 0.1 w1 + 0.1 w2) }
#' \deqn{ \alpha2 = exp(2 + 0.3 w2 + 0.3 w3) }
#' \deqn{ \alpha3 = exp(1.5 + 0.2 w1 + 0.1 w2 + 0.1 w3) }
#' \deqn{ \beta = exp(0.7 + 0.1 w1 + 0.1 w2 + 0.2 w3) }
#' For each set of parameters \eqn{(\alpha1, \alpha2, \alpha3, \beta)} a random sample from this bivariate gamma distribution is simulated, i.e. \eqn{y ~ BG(\alpha1, \alpha2, \alpha3, \beta)}. The gating is simulated from a logistic regression model where the regression coeffcients are
#' \deqn{logit(p) = 10 + 40 w1 + 30 w2 + 100 w3 ,}
#' based on which a binomial simulation is used, and the two components are sampled from the simulations above to form the data set. The final simulated data set consists of 232 samples from component 1 and 268 observations from component 2.
#'
#' @source Hu, S., Murphy, T. B. and O'Hagan, A. (2019) Bivariate Gamma Mixture of Experts Models for Joint Claims Modeling. To appear.
"fullsim"

#' @rdname fullsim
"fullsim.test"
