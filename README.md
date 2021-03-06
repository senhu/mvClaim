
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mvClaim R package

Multivariate general insurance claims modelling

Authors:

  - Sen Hu (<ethansen.hu@gmail.com>)
  - Thomas Brendan Murphy (<brendan.murphy@ucd.ie>)
  - Adrian O’Hagan (<adrian.ohagan@ucd.ie>)

## Description

The `mvClaim` package provides a flexible modelling framework of mixture
of experts (MoE) using bivariate gamma distributions, as introduced in
Hu et al (2019). It implements the bivariate gamma distribution proposed
by Cheriyan (1941) and Ramabhadran (1951), which has not received much
attention in the past. Depends on the parameterization of the framework,
its interpretation, and model parsimony issue, a family of models are
implemented including:

  - Bivariate gamma distribution estimation
  - Bivariate gamma regression
  - Mixture of bivariate gamma clustering
  - Mixture of bivariate gamma regressions (i.e. mixture of bivariate
    gamma clustering with covariates)

The `mvClaim` package also provides a similar mixture model framework of
finite mixture of copula regressions with gamma as marginal
distribution, including

  - Copula regression with gamma marginal distributions
  - Mixture of copula regressions with gamma marginal distributions
    (i.e. mixture of copulas clustering with covariates)

## Installation

You can install the latest development version of `mvClaim` from GitHub:

``` r
install.packages("devtools")
devtools::install_github("senhu/mvClaim")
```

Then the package can be loaded with:

``` r
library(mvClaim)
```

## Bivariate Gamma MoE Example

This README file follows a package vignette format, and an example is
briefly demonstrated using a simulated data set called `gatingsim` as
illustrated in Hu et al (2019). The data were simulated based on a
gating network MoE, i.e. covariates were used only in the gating network
to assist with identifying which component the observation was from.

``` r
data("gatingsim")
```

First without considering covariates, assuming there is only one
bivariate gamma distribution and no mixtures, the distribution
parameters can be estimated via

``` r
dist.est <- BGE(gatingsim[,1:2])
```

Still without considering covariates, if assuming that the data are a
mixture of two bivariate gamma distributions, there are three model
types within this bivariate gamma MoE family: “CC”, “CI” and “IC”.
Suppose the model type is “CC”, the mixture model can be fitted via:

``` r
mod1 <- MBGC(modelName = "CC", y=c("y1","y2"),
             G=2, gating = "C", data=gatingsim)
```

When taking covariates into consideration, we know the true data
generating process for the data set `gatingsim`, based on which the
optimal model can be fitted via

``` r
mod2 <- MBGC(modelName = "CC", y=c("y1","y2"),
             G=2, gating = ~w1+w2+w3, data=gatingsim)
```

In cases when covariates enter all gating and expert networks such as
the case of the artificial data set `fullsim`, a model can be fitted via

``` r
mod3 <- MBGR(modelName = "VV", y=c("y1","y2"),
             G=2, data = fullsim,
             f1     = ~ w1 + w2 + w3,
             f2     = ~ w1 + w2 + w3,
             f3     = ~ w1 + w2 + w3,
             f4     = ~ w1 + w2 + w3,
             gating = ~ w1 + w2 + w3)
```

In such cases, there are 9 different model types: “VC”, “VI”, “VV”,
“VE”, “CV”, “IV”, “EV”, “EC”, “CE”. Typically model selection issue
needs to be addressed if the true model is unknown. We recommend a
stepwise forward selection starting from a null model, i.e. fitting one
bivariate gamma distribution without covariates. More details can be
found in Hu et al (2019).

## Mixture of Copula Regressions Example

Another example is briefly demonstrated here using a simulated data set
called `simdat.mcgr` as illustrated in Hu and O’Hagan (2019). The data
were simulated based on a mixture of two copula regressions,
i.e. covariates were used to estimate the marginal distributions to
assist the mixture of copulas estimation. Details of the sample data
simulation process can be found in Hu and O’Hagan (2019).

``` r
data("simdat.mcgr")
```

Because we know the true data generating process for the data set
`simdat.mcgr`, based on which the model can be fitted via

``` r
mod4 <- MCGR(copula = list(gumbelCopula(dim=2), frankCopula(dim=2)),
             f1 = y1 ~ x1+x2,
             f2 = y2 ~ x1+x2,
             G = 2,
             gating = "C"
             data=simdat.mcgr)
```

If `G=1`, as an alternative, copula regression with gamma marginals can
be fitted via, for example:

``` r
mod5 <- copreg.gamma(f1 = y1 ~ x1+x2,
                     f2 = y2 ~ x1+x2,
                     copula = gumbelCopula(dim=2),
                     data = simdat.mcgr)
```

## Reference

1.  Hu, S., Murphy, T. B. and O’Hagan, A. (2019) Bivariate Gamma Mixture
    of Experts Models for Joint Claims Modeling. To appear.

2.  Hu, S. and O’Hagan, A. (2019) Copula Averaging for Tail Dependence
    in General Insurance Claims Data. To appear.
