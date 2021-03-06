#' Plot Mixture of Bivariate Gamma Regressions and Clustering Results
#'
#' Plot mixture of bivariate gamma regressions, or model-based clustering with bivariage gamma distributions and covariates: classification, uncertainty
#'
#' @param x Output from MBGR functions.
#' @param what The type of graph requested, either \code{"classification"}, or \code{"uncertainty"}, or \code{"fitted"}.
#' @param xlab,ylab Optional labels for the x-axis and the y-axis.
#' @param col The colors of points in clusters.
#' @param pch The plotting characters or symbols in clusters.
#' @param main A logical or \code{NULL} indicating whether or not to add a title to the plot identifying the type of plot.
#' @param ... Other graphical parameters.
#'
#' @return For BGR model, fitted values are plotted.
#'
#'   For MBGR model:
#'   \item{what=="classification"}{a plot showing the clustering labels.}
#'   \item{what=="uncertainty"}{a plot of classification uncertainty.}
#'   \item{what=="fitted"}{a plot of fitted values.}
#'   For MBGC model:
#'   \item{what=="classification"}{a plot showing the clustering labels.}
#'   \item{what=="uncertainty"}{a plot of classification uncertainty.}
#'
#' @examples
#' \donttest{
#' m1 <- MBGR(modelName = "CE", y=c("y1","y2"),
#'            data   = fullsim, G=2,
#'            f1     = ~ w1 + w2,
#'            f2     = ~ w2 + w3,
#'            f3     = ~ w1 + w2 + w3,
#'            f4     = ~ w1 + w2 + w3,
#'            gating = "C", verbose = FALSE)
#' plot.MBGR(m1, what="classification")
#' plot.MBGR(m1, what="uncertainty")
#' plot.MBGR(m1, what="fitted")
#' m2 <- BGR(modelName = "EI",
#'           y = c("y1","y2"), data = fullsim,
#'           f1     = ~ w1 + w2,
#'           f2     = ~ w2 + w3,
#'           f3     = ~ w1 + w2 + w3, verbose = FALSE)
#' plot.BGR(m2)
#' }
#'
#' @importFrom graphics points plot
#' @export plot.MBGR
#' @export

plot.MBGR <- function(x, what="classification",
                      col=NULL,
                      pch=NULL,
                      xlab=NULL,
                      ylab=NULL,
                      main=FALSE, ...){
  object <- x
  if (!inherits(object, "MBGR"))
    stop("object not of class \"MBGR\"")
  data <- object$y
  n    <- object$n
  G    <- object$G
  z    <- object$z
  class<- object$class
  newxlab <- if (is.null(xlab)) colnames(data)[1] else xlab
  newylab <- if (is.null(ylab)) colnames(data)[2] else ylab
  main.title <- if ( isTRUE(main) ) what else if (is.character(main)) main else NULL
  if (is.null(col)) points.col <- object$class else {
    points.col <- rep(NULL, n)
    for (gg in 1:G){ points.col <- replace(points.col,(object$class==gg),col[gg]) }
  }
  switch (what,
          classification = {
            #points.options <- c(0,16,3,2,23,4,11,5,15,1,6,7,8,9,10,12,13,14,17,18,19,20,21,22,24,25)
            #points.pch <- if (is.null(pch)) points.options[object$class] else pch
            if (is.null(pch)) points.pch <- c(0:20)[object$class] else {
              points.pch <- rep(NULL, n)
              for (gg in 1:G){ points.pch <- replace(points.pch,(object$class==gg),pch[gg]) }
            }
            plot(data, col=points.col, pch=points.pch,
                 main=main.title, xlab=newxlab, ylab=newylab, ...)
          },
          uncertainty = {
            plot(data,type="n",main=main.title, xlab=newxlab, ylab=newylab)
            for (gg in 1:G) points(data[class==gg,], col=points.col[class==gg],cex=3*(1-z[class==gg,gg]), pch=16)
          },
          fitted = {
            plot(rbind(data,object$fitted.values),
                 type="n",main=main.title,xlab=newxlab,ylab=newylab)
            points(data,col=points.col,pch=20,cex=.6)
            points(object$fitted.values, pch=15, col="brown")
          },
          stop("invalid plot type.")
  )
}

#' @rdname plot.MBGR
#' @export plot.BGR
#' @export

plot.BGR <- function(x,
                     col=NULL,
                     pch=NULL,
                     xlab=NULL,
                     ylab=NULL,
                     main=FALSE, ...){
  object <- x
  if (!inherits(object, "BGR"))
    stop("object not of class \"BGR\"")
  data <- object$y
  n    <- object$n
  newxlab <- if (is.null(xlab)) colnames(data)[1] else xlab
  newylab <- if (is.null(ylab)) colnames(data)[2] else ylab
  main.title <- if ( isTRUE(main) ) "BGR fitted values" else if (is.character(main)) main else NULL
  plot(rbind(data,object$fitted.values), type="n",
       main=main.title,xlab=newxlab,ylab=newylab,...)
  points(data,pch=20,cex=.6)
  points(object$fitted.values, pch=15, col="brown")
}


#' @rdname plot.MBGR
#' @export plot.MBGC
#' @export

plot.MBGC <- function(x, what="classification",
                      col=NULL,
                      pch=NULL,
                      xlab=NULL,
                      ylab=NULL,
                      main=FALSE, ...){
  object <- x
  if (!inherits(object, "MBGC"))
    stop("object not of class \"MBGC\"")
  data <- object$y
  n    <- object$n
  G    <- object$G
  z    <- object$z
  class<- object$class
  newxlab <- if (is.null(xlab)) colnames(data)[1] else xlab
  newylab <- if (is.null(ylab)) colnames(data)[2] else ylab
  main.title <- if ( isTRUE(main) ) what else if (is.character(main)) main else NULL
  if (is.null(col)) points.col <- object$class else {
    points.col <- rep(NULL, n)
    for (gg in 1:G){ points.col <- replace(points.col,(object$class==gg),col[gg]) }
  }
  switch (what,
          classification = {
            #points.options <- c(0,16,3,2,23,4,11,5,15,1,6,7,8,9,10,12,13,14,17,18,19,20,21,22,24,25)
            #points.pch <- if (is.null(pch)) c(0:20)[object$class] else pch
            if (is.null(pch)) points.pch <- c(0:20)[object$class] else {
              points.pch <- rep(NULL, n)
              for (gg in 1:G){ points.pch <- replace(points.pch,(object$class==gg),pch[gg]) }
            }
            plot(data, col=points.col, pch=points.pch,
                 main=main.title, xlab=newxlab, ylab=newylab, ...)
          },
          uncertainty = {
            plot(data,type="n",main=main.title, xlab=newxlab, ylab=newylab)
            for (gg in 1:G) points(data[class==gg,], col=points.col[class==gg],cex=3*(1-z[class==gg,gg]), pch=16)
          },
          stop("invalid plot type.")
  )
}
