BGR <- function(modelName = c("EE","EI","IE"),
                response,
                data,
                l1,
                l2,
                l3,
                l4,
                expo=NULL,
                maxit=10,
                tol=1e-5,
                Aitken = TRUE,
                verbose=TRUE){
  if (modelName == "EE"){
    # BGR.EE: all alpha's and beta are regressed on covariates
    if (is.null(l1)){ stop("l1 must be supplied for EE model type") }
    if (is.null(l2)){ stop("l2 must be supplied for EE model type") }
    if (is.null(l3)){ stop("l3 must be supplied for EE model type") }
    if (is.null(l4)){ stop("l4 must be supplied for EE model type") }
    return(BGR.EE(data=data, response=response,
                  l1=l1, l2=l2, l3=l3, l4=l4,
                  expo=expo, maxit=maxit, tol=tol,
                  Aitken = Aitken, verbose=verbose))
  }
  if (modelName == "EI"){
    # BGR.EI: same beta for all observations
    if (is.null(l1)){ stop("l1 must be supplied for EI model type") }
    if (is.null(l2)){ stop("l2 must be supplied for EI model type") }
    if (is.null(l3)){ stop("l3 must be supplied for EI model type") }
    return(BGR.EI(data=data, response=response,
                  l1=l1, l2=l2, l3=l3,
                  expo=expo, maxit=maxit, tol=tol,
                  Aitken=Aitken, verbose=verbose))
  }
  if (modelName == "IE"){
    # BGR.IE: only beta is regressed on its covariates, not alpha
    if (is.null(l4)){ stop("l4 must be supplied for IE model type") }
    return(BGR.IE(data=data, response=response,
                  l4=l4,
                  expo=expo, maxit=maxit, tol=tol,
                  Aitken=Aitken, verbose=verbose))
  }
}









