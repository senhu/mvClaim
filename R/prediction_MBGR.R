predict.mbgr4 <- function(object, newdata){
  
  l1.n    <- object$formula[[1]]
  matrix1 <- as.matrix(Matrix::sparse.model.matrix(l1.n, data=newdata))
  l2.n    <- object$formula[[2]]
  matrix2 <- as.matrix(Matrix::sparse.model.matrix(l2.n, data=newdata))
  l3.n    <- object$formula[[3]]
  matrix3 <- as.matrix(Matrix::sparse.model.matrix(l3.n, data=newdata))
  l4.n    <- object$formula[[4]]
  matrix4 <- as.matrix(Matrix::sparse.model.matrix(l4.n, data=newdata))
  lp.n    <- object$formula[[5]]

  m      <- dim(newdata)[1]
  
  a1 <- exp(matrix1 %*% object$coefficients[[1]])
  a2 <- exp(matrix2 %*% object$coefficients[[2]])
  a3 <- exp(matrix3 %*% object$coefficients[[3]])
  b  <- exp(matrix4 %*% object$coefficients[[4]])
  
  if (is.null(lp.n)){
    pz         <- object$clusterProb
    pred1.test <- ((a1 + a3) / b) %*% pz
    pred2.test <- ((a2 + a3) / b) %*% pz
  } else {
    matrixp    <- as.matrix(Matrix::sparse.model.matrix(lp.n, data=newdata))
    denom      <- 1 + rowSums(exp(matrixp %*% t(object$coefficients[[5]])))
    numerator  <- cbind(rep(1, m), exp(matrixp %*% t(object$coefficients[[5]])) )
    pz         <- sweep(numerator, 1, denom, FUN="/")
    pred1.test <- rowSums(pz * ((a1 + a3)/b) )
    pred2.test <- rowSums(pz * ((a2 + a3)/b) )
  }
  
  pred.res <- cbind(pred1.test, pred2.test)
  colnames(pred.res) <- NULL
  
  result <- list(fit        = pred.res,
                 alpha1.fit = a1,
                 alpha2.fit = a2,
                 alpha3.fit = a3, 
                 beta.fit   = b,
                 tau.fit    = pz)
  return(result)
}



predict.mbgr2 <- function(object, newdata){
  
  l1.n    <- object$formula[[1]]
  matrix1 <- as.matrix(Matrix::sparse.model.matrix(l1.n, data=newdata))
  l2.n    <- object$formula[[2]]
  matrix2 <- as.matrix(Matrix::sparse.model.matrix(l2.n, data=newdata))
  l3.n    <- object$formula[[3]]
  matrix3 <- as.matrix(Matrix::sparse.model.matrix(l3.n, data=newdata))
  lp.n       <- object$formula[[4]]
  
  m      <- dim(newdata)[1]
  
  a1 <- exp(matrix1 %*%  object$coefficients[[1]])
  a2 <- exp(matrix2 %*%  object$coefficients[[2]])
  a3 <- exp(matrix3 %*%  object$coefficients[[3]])
  b  <- object$beta
  
  if (is.null(lp.n)){
    pz         <- object$clusterProb
    pred1.test <- sweep((a1+a3),2,b,FUN="/") %*% pz
    pred2.test <- sweep((a2+a3),2,b,FUN="/") %*% pz
  } else {
    matrixp    <- as.matrix(Matrix::sparse.model.matrix(lp.n, data=newdata))
    
    denom      <- 1 + rowSums(exp(matrixp %*% t(object$coefficients[[4]])))
    numerator  <- cbind(rep(1, m), exp(matrixp %*% t(object$coefficients[[4]])) )
    pz         <- sweep(numerator, 1, denom, FUN="/")
    pred1.test <- rowSums(pz * sweep((a1+a3),2,b,FUN="/") )
    pred2.test <- rowSums(pz * sweep((a2+a3),2,b,FUN="/") )
  }
  
  pred.res <- cbind(pred1.test, pred2.test)
  colnames(pred.res) <- NULL
  
  result <- list(fit        = pred.res,
                 alpha1.fit = a1,
                 alpha2.fit = a2,
                 alpha3.fit = a3, 
                 beta.fit   = b,
                 tau.fit    = pz)
  
  return(result)
}






