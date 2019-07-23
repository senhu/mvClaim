
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mvClaim R package

Mixture of experts model with bivariate gamma distributions

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

## Artificial Data Example

This README file follows a package vignette format, and an example is
briefly demonstrated using a simulated data set called `gatingsim` as
illustrated in Hu et al (2019). The data were simulated based on a
gating network MoE, i.e. covariates were used only in the gating network
to assist with identifying which component the observation was from. The
pairwise plot of the data is shown below.
<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />

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
generating process, based on which the optimal model can be fitted via

``` r
mod2 <- MBGC(modelName = "CC", y=c("y1","y2"),
             G=2, gating = ~w1+w2+w3, data=gatingsim)
```

Although it is not the case that for the true data generating process
covariates only appear in the gating network, in cases when covariates
also enter all gating and expert networks, a model can be fitted via

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

## Reference

Hu, S., Murphy, T. B. and O’Hagan, A. (2019) Bivariate Gamma Mixture of
Experts Models for Joint Claims Modeling. To appear.
