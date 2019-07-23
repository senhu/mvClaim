expected.latent <- function(y, alpha, beta){
  y1 <- y[1]
  y2 <- y[2]
  alpha1 <- alpha[1]
  alpha2 <- alpha[2]
  alpha3 <- alpha[3]
  beta <- beta

  integrand.numerator.logs <- function(x3){
    log(x3) * dgamma(y1-x3,shape=alpha1,rate=beta, log=FALSE)*
      dgamma(y2-x3,shape=alpha2,rate=beta, log=FALSE)*
      dgamma(x3,shape=alpha3,rate=beta, log=FALSE)
  }


  integrand.numerator.logy1s <- function(x3){
    log(y1-x3) * dgamma(y1-x3,shape=alpha1,rate=beta, log=FALSE)*
      dgamma(y2-x3,shape=alpha2,rate=beta, log=FALSE)*
      dgamma(x3,shape=alpha3,rate=beta, log=FALSE)
  }

  integrand.numerator.logy2s <- function(x3){
    log(y2-x3) * dgamma(y1-x3,shape=alpha1,rate=beta, log=FALSE)*
      dgamma(y2-x3,shape=alpha2,rate=beta, log=FALSE)*
      dgamma(x3,shape=alpha3,rate=beta, log=FALSE)
  }

  #-------------------------
  int.numerator.s.res <- log(alpha3)-log(beta)+dbivgamma(y, c(alpha1, alpha2, alpha3+1), beta, log = TRUE)$value
  #-------------------------
  int.numerator.logs <- distrEx::distrExIntegrate(integrand.numerator.logs,
                                                  lower=0, upper=min(y1,y2),
                                                  rel.tol = 1e-6, order=1e+4)

  if (abs(int.numerator.logs)!=Inf){int.numerator.logs.res <- int.numerator.logs }
  if (int.numerator.logs==Inf){
    int.numerator.logs <- distrEx::distrExIntegrate(integrand.numerator.logs,
                                                    lower=.Machine$double.eps^{1/3},
                                                    upper=(min(y1, y2)-.Machine$double.eps^{1/3}),
                                                    rel.tol = 1e-6,
                                                    order = 1e+4)
    if (is.na(int.numerator.logs) || is.nan(int.numerator.logs) || abs(int.numerator.logs)==Inf ){
      print(int.numerator.logs); stop("error in expected.latent: numerator.logs")
    } else int.numerator.logs.res <- int.numerator.logs
  }
  #-------------------------
  int.numerator.logy1s <- distrEx::distrExIntegrate(integrand.numerator.logy1s,
                                                    lower=0, upper=min(y1,y2),
                                                    rel.tol = 1e-6, order=1e+4)
  if (abs(int.numerator.logy1s)!=Inf){int.numerator.logy1s.res <- int.numerator.logy1s}
  if (abs(int.numerator.logy1s)==Inf){
    int.numerator.logy1s <- distrEx::distrExIntegrate(integrand.numerator.logy1s,
                                                      lower=.Machine$double.eps^{1/3},
                                                      upper=(min(y1, y2)-.Machine$double.eps^{1/3}),
                                                      rel.tol = 1e-6,
                                                      order = 1e+5)
    if (is.na(int.numerator.logy1s) || is.nan(int.numerator.logy1s) || abs(int.numerator.logy1s)==Inf){
      print(int.numerator.logy1s); stop("error in expected.latent: numerator.logy1s")
    } else int.numerator.logy1s.res <- int.numerator.logy1s
  }
  #---------------------------
  int.numerator.logy2s <- distrEx::distrExIntegrate(integrand.numerator.logy2s,
                                                    lower=0, upper=min(y1,y2),
                                                    rel.tol = 1e-6, order=1e+4)
  if (abs(int.numerator.logy2s)!=Inf){int.numerator.logy2s.res <- int.numerator.logy2s }
  if (abs(int.numerator.logy2s)==Inf){
    int.numerator.logy2s <- distrEx::distrExIntegrate(integrand.numerator.logy2s,
                                                      lower=.Machine$double.eps^{1/3},
                                                      upper=(min(y1, y2)-.Machine$double.eps^{1/3}),
                                                      rel.tol = 1e-6,
                                                      order = 1e+5)
    if (is.na(int.numerator.logy2s) || is.nan(int.numerator.logy2s) || abs(int.numerator.logy2s) == Inf){
      print(int.numerator.logy2s); stop("error in expected.latent: numerator.logy2s")
    } else int.numerator.logy2s.res <- int.numerator.logy2s
  }
  #---------------------------
  int.denominator.res <- dbivgamma(y, alpha, beta, log=TRUE)$value
  #------------------------------

  Expected.s      <- exp(int.numerator.s.res - int.denominator.res)
  Expected.logs   <- int.numerator.logs.res/exp(int.denominator.res)
  Expected.logy1s <- int.numerator.logy1s.res/exp(int.denominator.res)
  Expected.logy2s <- int.numerator.logy2s.res/exp(int.denominator.res)

  if (exp(int.denominator.res) == 0){
    Expected.logs   <- log(Expected.s)
    Expected.logy1s <- log(y1-Expected.s)
    Expected.logy2s <- log(y2-Expected.s)
  }

  if (is.nan(Expected.s)){stop("NaN value returned for Expected.s")}
  if (is.na(Expected.s)){stop("NA value returned for Expected.s")}
  if (is.nan(Expected.logs)){stop("NaN value returned for Expected.logs")}
  if (is.na(Expected.logs)){stop("NA value returned for Expected.logs")}
  if (is.nan(Expected.logy1s)){stop(cat("NaN value returned for Expected.logy1s",
                                        "numerator=", int.numerator.logy1s.res,
                                        "  denominator = exp of ", int.denominator.res,"\n"))}
  if (is.na(Expected.logy1s)){stop("NA value returned for Expected.logy1s")}
  if (is.nan(Expected.logy2s)){stop(cat("NaN value returned for Expected.logy2s",
                                        "numerator=", int.numerator.logy2s.res,
                                        "  denominator = exp of ", int.denominator.res, "\n"))}
  if (is.na(Expected.logy2s)){stop("NA value returned for Expected.logy2s")}

  return(list("Expected.s"=Expected.s,
              "Expected.logs"=Expected.logs,
              "Expected.logy1s"=Expected.logy1s,
              "Expected.logy2s"=Expected.logy2s,
              "s.numerator"=int.numerator.s.res,
              "logs.numerator"=int.numerator.logs,
              "logy1s.numerator"=int.numerator.logy1s,
              "logy2s.numerator"=int.numerator.logy2s,
              "denominator"=int.denominator.res))
}

