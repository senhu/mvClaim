#' Predict Method for Bivariate Gamma Regressions
#'
#' Obtains predictions from the fitted bivariate gamma regression models.
#'
#' @param object a fitted model object of class \code{"BGR"} or \code{"MBGR"}
#' @param newdata a data frame that prediction is based upon
#' @param ... further arguments passed to or from other methods.
#'
#' @return a list with components (depends on model class):
#'   \item{fit}{predictions.}
#'   \item{alpha1.fit}{predicted alpha1 values.}
#'   \item{alpha2.fit}{predicted alpha2 values.}
#'   \item{alpha3.fit}{predicted alpha3 values.}
#'   \item{beta.fit}{predicted beta values.}
#'   \item{tau.fit}{predicted mixing proportion values, if covariates enter the gating network.}
#'
#' @examples
#' mod1 <- BGR(modelName = "EE",
#'             y = c("y1","y2"), data = fullsim,
#'             f1 = ~ w1 + w2,
#'             f2 = ~ w2 + w3,
#'             f3 = ~ w1 + w2 + w3,
#'             f4 = ~ w1 + w2 + w3,
#'             verbose= FALSE)
#' fitted1 <- predict(mod1, newdata=fullsim)
#' mod2 <- BGR(modelName = "EI",
#'             y = c("y1","y2"), data = fullsim,
#'             f1     = ~ w1 + w2,
#'             f2     = ~ w2 + w3,
#'             f3     = ~ w1 + w2 + w3,
#'             verbose= FALSE)
#' fitted2 <- predict(mod2, newdata=fullsim)
#' mod3 <- BGR(modelName = "IE",
#'             y = c("y1","y2"), data = fullsim,
#'             f4     = ~ w1 + w2 + w3,
#'             verbose= FALSE)
#' fitted3 <- predict(mod3, newdata=fullsim)
#' plot(fitted1$fit)
#' plot(fitted2$fit)
#' plot(fitted3$fit)

predict.BGR <- function(object, newdata, ...){
  modelName <- object$modelName
  switch(modelName,
         EE = {    # BGR.EE: all alpha's and beta are regressed on covariates
           return(predict.BGR.EE(object=object, newdata=newdata, ...))
         },
         EI = {    # BGR.EI: same beta for all observations
           return(predict.BGR.EI(object=object, newdata=newdata, ...))
         },
         IE = {# BGR.IE: only beta is regressed on its covariates, not alpha
           return(predict.BGR.IE(object=object, newdata=newdata, ...))
         },
         stop("invalid model type name")
  )
}

predict.BGR.EE <- function(object, newdata, ...){

  l1.n    <- object$formula[[1]]
  matrix1 <- as.matrix(stats::model.matrix(l1.n, data=newdata))
  l2.n    <- object$formula[[2]]
  matrix2 <- as.matrix(stats::model.matrix(l2.n, data=newdata))
  l3.n    <- object$formula[[3]]
  matrix3 <- as.matrix(stats::model.matrix(l3.n, data=newdata))
  l4.n    <- object$formula[[4]]
  matrix4 <- as.matrix(stats::model.matrix(l4.n, data=newdata))

  a1 <- exp(matrix1 %*% object$coefficients[[1]])
  a2 <- exp(matrix2 %*% object$coefficients[[2]])
  a3 <- exp(matrix3 %*% object$coefficients[[3]])
  b  <- exp(matrix4 %*% object$coefficients[[4]])

  pred1.test <- rowSums((a1 + a3)/b)
  pred2.test <- rowSums((a2 + a3)/b)
  pred.res <- cbind(pred1.test, pred2.test)
  colnames(pred.res) <- NULL
  result <- list(fit        = pred.res,
                 alpha1.fit = a1,
                 alpha2.fit = a2,
                 alpha3.fit = a3,
                 beta.fit   = b)
  return(result)
}

predict.BGR.EI <- function(object, newdata, ...){

  l1.n    <- object$formula[[1]]
  matrix1 <- as.matrix(stats::model.matrix(l1.n, data=newdata))
  l2.n    <- object$formula[[2]]
  matrix2 <- as.matrix(stats::model.matrix(l2.n, data=newdata))
  l3.n    <- object$formula[[3]]
  matrix3 <- as.matrix(stats::model.matrix(l3.n, data=newdata))

  a1 <- exp(matrix1 %*% object$coefficients[[1]])
  a2 <- exp(matrix2 %*% object$coefficients[[2]])
  a3 <- exp(matrix3 %*% object$coefficients[[3]])
  b  <- object$beta

  pred1.test <- rowSums((a1 + a3)/b)
  pred2.test <- rowSums((a2 + a3)/b)
  pred.res <- cbind(pred1.test, pred2.test)
  colnames(pred.res) <- NULL
  result <- list(fit        = pred.res,
                 alpha1.fit = a1,
                 alpha2.fit = a2,
                 alpha3.fit = a3,
                 beta.fit   = b)
  return(result)
}

predict.BGR.IE <- function(object, newdata, ...){
  l4.n    <- object$formula
  matrix4 <- as.matrix(stats::model.matrix(l4.n, data=newdata))
  a1 <- object$alpha1
  a2 <- object$alpha2
  a3 <- object$alpha3
  b  <- exp(matrix4 %*% object$coefficients)
  pred1.test <- rowSums((a1 + a3)/b)
  pred2.test <- rowSums((a2 + a3)/b)
  pred.res <- cbind(pred1.test, pred2.test)
  colnames(pred.res) <- NULL
  result <- list(fit        = pred.res,
                 alpha1.fit = a1,
                 alpha2.fit = a2,
                 alpha3.fit = a3,
                 beta.fit   = b)
  return(result)
}






