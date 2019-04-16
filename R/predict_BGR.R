

predict <- function(object, newdata){
  UseMethod("predict", object)
}

predict.BGR.EE <- function(object,
                           newdata){

  l1.n    <- object$formula[[1]]
  matrix1 <- as.matrix(Matrix::sparse.model.matrix(l1.n, data=newdata))
  l2.n    <- object$formula[[2]]
  matrix2 <- as.matrix(Matrix::sparse.model.matrix(l2.n, data=newdata))
  l3.n    <- object$formula[[3]]
  matrix3 <- as.matrix(Matrix::sparse.model.matrix(l3.n, data=newdata))
  l4.n    <- object$formula[[4]]
  matrix4 <- as.matrix(Matrix::sparse.model.matrix(l4.n, data=newdata))

  a1 <- exp(matrix1 %*% object$coef[[1]])
  a2 <- exp(matrix2 %*% object$coef[[2]])
  a3 <- exp(matrix3 %*% object$coef[[3]])
  b  <- exp(matrix4 %*% object$coef[[4]])

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



predict.BGR.EI <- function(object,
                           newdata){

  l1.n    <- object$formula[[1]]
  matrix1 <- as.matrix(Matrix::sparse.model.matrix(l1.n, data=newdata))
  l2.n    <- object$formula[[2]]
  matrix2 <- as.matrix(Matrix::sparse.model.matrix(l2.n, data=newdata))
  l3.n    <- object$formula[[3]]
  matrix3 <- as.matrix(Matrix::sparse.model.matrix(l3.n, data=newdata))

  a1 <- exp(matrix1 %*% object$coef[[1]])
  a2 <- exp(matrix2 %*% object$coef[[2]])
  a3 <- exp(matrix3 %*% object$coef[[3]])
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

predict.BGR.IE <- function(object,
                           newdata){

  l4.n    <- object$formula[[1]]
  matrix4 <- as.matrix(Matrix::sparse.model.matrix(l4.n, data=newdata))

  a1 <- object$alpha1
  a2 <- object$alpha2
  a3 <- object$alpha3
  b  <- exp(matrix4 %*% object$coef[[4]])

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
