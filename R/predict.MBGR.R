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
#' m1 <- MBGR(modelName = "VC", G=2,
#'            y=c("y1","y2"), data = fullsim,
#'            f1     = ~ w1 + w2,
#'            f2     = ~ w2 + w3,
#'            f3     = ~ w1 + w2 + w3,
#'            f4     = ~ w1 + w2 + w3,
#'            gating = "C", verbose=FALSE)
#' fitted1 <- predict(m1, newdata=fullsim)


predict.MBGR <- function(object, newdata, ...){
  modelName <- object$modelName
  switch(modelName,
         VC = {
           res <- predict.MBGR.VC(object=object, newdata=newdata, ...)
         },
         CV = {
           res <- predict.MBGR.CV(object=object, newdata=newdata, ...)
         },
         VV = {
           res <- predict.MBGR.VV(object=object, newdata=newdata, ...)
         },
         VI = {
           res <- predict.MBGR.VI(object=object, newdata=newdata, ...)
         },
         IV = {
           res <- predict.MBGR.IV(object=object, newdata=newdata, ...)
         },
         VE = {
           res <- predict.MBGR.VE(object=object, newdata=newdata, ...)
         },
         EV = {
           res <- predict.MBGR.EV(object=object, newdata=newdata, ...)
         },
         EC = {
           res <- predict.MBGR.EC(object=object, newdata=newdata, ...)
         },
         CE = {
           res <- predict.MBGR.CE(object=object, newdata=newdata, ...)
         } )
  return(res)
}



# MBGR2
predict.MBGR.VC <- function(object, newdata){
  l1.n    <- object$formula[[1]]
  matrix1 <- as.matrix(stats::model.matrix(l1.n, data=newdata))
  l2.n    <- object$formula[[2]]
  matrix2 <- as.matrix(stats::model.matrix(l2.n, data=newdata))
  l3.n    <- object$formula[[3]]
  matrix3 <- as.matrix(stats::model.matrix(l3.n, data=newdata))
  lp.n    <- object$formula[[4]]
  m      <- dim(newdata)[1]
  a1 <- exp(matrix1 %*%  object$coefficients$expert[[1]])
  a2 <- exp(matrix2 %*%  object$coefficients$expert[[2]])
  a3 <- exp(matrix3 %*%  object$coefficients$expert[[3]])
  b  <- object$beta
  if (class(lp.n)!="formula"){
    pz         <- object$pro
    pred1.test <- sweep((a1+a3),2,b,FUN="/") %*% pz
    pred2.test <- sweep((a2+a3),2,b,FUN="/") %*% pz
  } else {
    matrixp    <- as.matrix(stats::model.matrix(lp.n, data=newdata))
    denom      <- 1 + rowSums(exp(matrixp %*% t(object$coefficients$gating)))
    numerator  <- cbind(rep(1, m), exp(matrixp %*% t(object$coefficients$gating)))
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
# MBGR3
predict.MBGR.CV <- function(object, newdata){
  l4.n    <- object$formula[[1]]
  matrix4 <- as.matrix(stats::model.matrix(l4.n, data=newdata))
  lp.n    <- object$formula[[2]]
  m      <- nrow(newdata)
  a1 <- object$alpha1
  a2 <- object$alpha2
  a3 <- object$alpha3
  b  <- exp(matrix4 %*% object$coefficients$expert[[1]])
  if (class(lp.n)!="formula"){
    pz         <- object$pro
    pred1.test <- sweep((1/b), 2, (a1+a3), FUN="*") %*% pz
    pred2.test <- sweep((1/b), 2, (a2+a3), FUN="*") %*% pz
  } else {
    matrixp    <- as.matrix(stats::model.matrix(lp.n, data=newdata))
    denom      <- 1 + rowSums(exp(matrixp %*% t(object$coefficients$gating)))
    numerator  <- cbind(rep(1, m), exp(matrixp %*% t(object$coefficients$gating)) )
    pz         <- sweep(numerator, 1, denom, FUN="/")
    pred1.test <- rowSums(pz * sweep((1/b), 2, (a1+a3), FUN="*"))
    pred2.test <- rowSums(pz * sweep((1/b), 2, (a2+a3), FUN="*"))
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
# MBGR4
predict.MBGR.VV <- function(object, newdata){
  l1.n    <- object$formula[[1]]
  matrix1 <- as.matrix(stats::model.matrix(l1.n, data=newdata))
  l2.n    <- object$formula[[2]]
  matrix2 <- as.matrix(stats::model.matrix(l2.n, data=newdata))
  l3.n    <- object$formula[[3]]
  matrix3 <- as.matrix(stats::model.matrix(l3.n, data=newdata))
  l4.n    <- object$formula[[4]]
  matrix4 <- as.matrix(stats::model.matrix(l4.n, data=newdata))
  lp.n       <- object$formula[[5]]
  m      <- nrow(newdata)
  a1 <- exp(matrix1 %*% object$coefficients$expert[[1]])
  a2 <- exp(matrix2 %*% object$coefficients$expert[[2]])
  a3 <- exp(matrix3 %*% object$coefficients$expert[[3]])
  b  <- exp(matrix4 %*% object$coefficients$expert[[4]])
  if (class(lp.n)!="formula"){
    pz         <- object$pro
    pred1.test <- ((a1 + a3) / b) %*% pz
    pred2.test <- ((a2 + a3) / b) %*% pz
  } else {
    matrixp    <- as.matrix(stats::model.matrix(lp.n, data=newdata))
    denom      <- 1 + rowSums(exp(matrixp %*% t(object$coefficients$gating)))
    numerator  <- cbind(rep(1, m), exp(matrixp %*% t(object$coefficients$gating)))
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
# MBGR5
predict.MBGR.VI <- function(object, newdata){
  l1.n    <- object$formula[[1]]
  matrix1 <- as.matrix(stats::model.matrix(l1.n, data=newdata))
  l2.n    <- object$formula[[2]]
  matrix2 <- as.matrix(stats::model.matrix(l2.n, data=newdata))
  l3.n    <- object$formula[[3]]
  matrix3 <- as.matrix(stats::model.matrix(l3.n, data=newdata))
  lp.n    <- object$formula[[4]]
  m       <- nrow(newdata)
  G       <- object$G
  a1 <- exp(matrix1 %*%  object$coefficients$expert[[1]])
  a2 <- exp(matrix2 %*%  object$coefficients$expert[[2]])
  a3 <- exp(matrix3 %*%  object$coefficients$expert[[3]])
  b  <- object$beta
  if (class(lp.n)!="formula"){
    pz         <- object$pro
    pred1.test <- sweep((a1+a3),2,rep(b,G),FUN="/") %*% pz
    pred2.test <- sweep((a2+a3),2,rep(b,G),FUN="/") %*% pz
  } else {
    matrixp    <- as.matrix(stats::model.matrix(lp.n, data=newdata))
    denom      <- 1 + rowSums(exp(matrixp %*% t(object$coefficients$gating)))
    numerator  <- cbind(rep(1, m), exp(matrixp %*% t(object$coefficients$gating)))
    pz         <- sweep(numerator, 1, denom, FUN="/")
    pred1.test <- rowSums(pz * sweep((a1+a3),2,rep(b,G),FUN="/") )
    pred2.test <- rowSums(pz * sweep((a2+a3),2,rep(b,G),FUN="/") )
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
# MBGR6
predict.MBGR.IV <- function(object, newdata){
  l4.n    <- object$formula[[1]]
  matrix4 <- as.matrix(stats::model.matrix(l4.n, data=newdata))
  lp.n    <- object$formula[[2]]
  m       <- nrow(newdata)
  G       <- object$G
  a1 <- object$alpha1
  a2 <- object$alpha2
  a3 <- object$alpha3
  b  <- exp(matrix4 %*% object$coefficients$expert[[1]])
  if (class(lp.n)!="formula"){
    pz         <- object$pro
    pred1.test <- sweep((1/b), 2, (rep(a1+a3,G)), FUN="*") %*% pz
    pred2.test <- sweep((1/b), 2, (rep(a2+a3,G)), FUN="*") %*% pz
  } else {
    matrixp    <- as.matrix(stats::model.matrix(lp.n, data=newdata))
    denom      <- 1 + rowSums(exp(matrixp %*% t(object$coefficients$gating)))
    numerator  <- cbind(rep(1, m), exp(matrixp %*% t(object$coefficients$gating)) )
    pz         <- sweep(numerator, 1, denom, FUN="/")
    pred1.test <- rowSums(pz * sweep((1/b), 2, (rep(a1+a3,G)), FUN="*"))
    pred2.test <- rowSums(pz * sweep((1/b), 2, (rep(a2+a3,G)), FUN="*"))
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
# MBGR7
predict.MBGR.VE <- function(object, newdata){
  l1.n    <- object$formula[[1]]
  matrix1 <- as.matrix(stats::model.matrix(l1.n, data=newdata))
  l2.n    <- object$formula[[2]]
  matrix2 <- as.matrix(stats::model.matrix(l2.n, data=newdata))
  l3.n    <- object$formula[[3]]
  matrix3 <- as.matrix(stats::model.matrix(l3.n, data=newdata))
  l4.n    <- object$formula[[4]]
  matrix4 <- as.matrix(stats::model.matrix(l4.n, data=newdata))
  lp.n    <- object$formula[[5]]
  m       <- dim(newdata)[1]
  G       <- object$G
  a1 <- exp(matrix1 %*% object$coefficients$expert[[1]])
  a2 <- exp(matrix2 %*% object$coefficients$expert[[2]])
  a3 <- exp(matrix3 %*% object$coefficients$expert[[3]])
  b  <- exp(matrix4 %*% object$coefficients$expert[[4]])
  if (class(lp.n)!="formula"){
    pz         <- object$pro
    pred1.test <- ((a1+a3) / rep.col(b,G)) %*% pz
    pred2.test <- ((a2+a3) / rep.col(b,G)) %*% pz
  } else {
    matrixp    <- as.matrix(stats::model.matrix(lp.n, data=newdata))
    denom      <- 1 + rowSums(exp(matrixp %*% t(object$coefficients$gating)))
    numerator  <- cbind(rep(1, m), exp(matrixp %*% t(object$coefficients$gating)))
    pz         <- sweep(numerator, 1, denom, FUN="/")
    pred1.test <- rowSums(pz * ((a1 + a3)/rep.col(b,G)))
    pred2.test <- rowSums(pz * ((a2 + a3)/rep.col(b,G)))
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
# MBGR8
predict.MBGR.EV <- function(object, newdata){
  l1.n    <- object$formula[[1]]
  matrix1 <- as.matrix(stats::model.matrix(l1.n, data=newdata))
  l2.n    <- object$formula[[2]]
  matrix2 <- as.matrix(stats::model.matrix(l2.n, data=newdata))
  l3.n    <- object$formula[[3]]
  matrix3 <- as.matrix(stats::model.matrix(l3.n, data=newdata))
  l4.n    <- object$formula[[4]]
  matrix4 <- as.matrix(stats::model.matrix(l4.n, data=newdata))
  lp.n    <- object$formula[[5]]
  m       <- dim(newdata)[1]
  G       <- object$G
  a1 <- exp(matrix1 %*% object$coefficients$expert[[1]])
  a2 <- exp(matrix2 %*% object$coefficients$expert[[2]])
  a3 <- exp(matrix3 %*% object$coefficients$expert[[3]])
  b  <- exp(matrix4 %*% object$coefficients$expert[[4]])
  if (class(lp.n)!="formula"){
    pz         <- object$pro
    pred1.test <- (rep.col(a1+a3,G) / b) %*% pz
    pred2.test <- (rep.col(a2+a3,G) / b) %*% pz
  } else {
    matrixp    <- as.matrix(stats::model.matrix(lp.n, data=newdata))
    denom      <- 1 + rowSums(exp(matrixp %*% t(object$coefficients$gating)))
    numerator  <- cbind(rep(1, m), exp(matrixp %*% t(object$coefficients$gating)))
    pz         <- sweep(numerator, 1, denom, FUN="/")
    pred1.test <- rowSums(pz * (rep.col(a1+a3,G)/b))
    pred2.test <- rowSums(pz * (rep.col(a2+a3,G)/b))
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
# MBGR9
predict.MBGR.EC <- function(object, newdata){
  l1.n    <- object$formula[[1]]
  matrix1 <- as.matrix(stats::model.matrix(l1.n, data=newdata))
  l2.n    <- object$formula[[2]]
  matrix2 <- as.matrix(stats::model.matrix(l2.n, data=newdata))
  l3.n    <- object$formula[[3]]
  matrix3 <- as.matrix(stats::model.matrix(l3.n, data=newdata))
  lp.n    <- object$formula[[4]]
  m       <- dim(newdata)[1]
  G       <- object$G
  a1 <- exp(matrix1 %*%  object$coefficients$expert[[1]])
  a2 <- exp(matrix2 %*%  object$coefficients$expert[[2]])
  a3 <- exp(matrix3 %*%  object$coefficients$expert[[3]])
  b  <- object$beta
  if (class(lp.n)!="formula"){
    pz         <- object$pro
    pred1.test <- sweep(rep.col(a1+a3,G),2,b,FUN="/") %*% pz
    pred2.test <- sweep(rep.col(a2+a3,G),2,b,FUN="/") %*% pz
  } else {
    matrixp    <- as.matrix(stats::model.matrix(lp.n, data=newdata))
    denom      <- 1 + rowSums(exp(matrixp %*% t(object$coefficients$gating)))
    numerator  <- cbind(rep(1, m), exp(matrixp %*% t(object$coefficients$gating)))
    pz         <- sweep(numerator, 1, denom, FUN="/")
    pred1.test <- rowSums(pz * sweep(rep.col(a1+a3,G),2,b,FUN="/"))
    pred2.test <- rowSums(pz * sweep(rep.col(a2+a3,G),2,b,FUN="/"))
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
#MBGR10
predict.MBGR.CE <- function(object, newdata){
  l4.n    <- object$formula[[1]]
  matrix4 <- as.matrix(stats::model.matrix(l4.n, data=newdata))
  lp.n    <- object$formula[[2]]
  m       <- nrow(newdata)
  G       <- object$G
  a1 <- object$alpha1
  a2 <- object$alpha2
  a3 <- object$alpha3
  b  <- exp(matrix4 %*% object$coefficients$expert[[1]])
  if (class(lp.n)!="formula"){
    pz         <- object$pro
    pred1.test <- sweep(rep.col(1/b,G),2,(a1+a3), FUN="*") %*% pz
    pred2.test <- sweep(rep.col(1/b,G),2,(a2+a3), FUN="*") %*% pz
  } else {
    matrixp    <- as.matrix(stats::model.matrix(lp.n, data=newdata))
    denom      <- 1 + rowSums(exp(matrixp %*% t(object$coefficients$gating)))
    numerator  <- cbind(rep(1, m), exp(matrixp %*% t(object$coefficients$gating)) )
    pz         <- sweep(numerator, 1, denom, FUN="/")
    pred1.test <- rowSums(pz * sweep(rep.col(1/b,G),2,(a1+a3), FUN="*"))
    pred2.test <- rowSums(pz * sweep(rep.col(1/b,G),2,(a2+a3), FUN="*"))
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


