
<!-- README.md is generated from README.Rmd. Please edit that file -->
mvClaim R package
=================

Mixture of Experts model with bivariate gamma and bivariate Poisson distributions

Description
-----------

The `mvClaim` package provides a flexible modelling framework of mixture of experts (MoE) using bivariate gamma distributions, as introduced in Hu et al. (2019). It implements the bivariate gamma distribution proposed by Cheriyan (1941) and Ramabhadran (1951), which has not received much attention in the past. Depends on the parameterization of the framework, its interpretation, and model parsimony issue, a family of models are implemented including:

-   Bivariate gamma distribution estimation
-   Bivariate gamma regression
-   Mixture of bivariate gamma clustering
-   Mixture of bivariate gamma regressions (i.e. mixture of bivariate gamma clustering with covariates)

Previous works in the literature have been investigated for bivariate and multivariate Poisson distribution which was constructed in a similar fashion, hence a framework of MoE using bivariate Poisson is also included in this package.

Installation
------------

You can install the latest development version of `mvClaim` from GitHub:

``` r
install.packages("devtools")
devtools::install_github("senhu/mvClaim")
```

Then the package can be loaded with:

``` r
library(mvClaim)
```

Artificial Data Example
-----------------------

Instead of a package vignette document, here an example is briefly demonstrated using a simulated data set called `gatingsim` as illustrated in the paper (Hu et al. 2019). The data were simulated based on a gating network MoE, i.e. covariates were used only in the gating network to assist with identifying which component the observation was from. The pairwise plot of the data is shown below. <img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />

First without considering covariates, assuming there is only one bivariate gamma distribution and no mixtures, the distribution parameters can be estimated by

``` r
dist.est <- BGE(gatingsim[,1:2])
```

or equivalently

``` r
dist.est <- MBGC(gatingsim[,1:2], G=1)
```

If assuming that the data are a mixture of two bivariate gamma distributions:

``` r
clust2 <- MBGC(gatingsim[,1:2], G=2)
```

When taking covariates into consideration, since we know the true data generating process, the optimal model can be fitted by

``` r
mod1 <- MBGR_Gating(response=c("y1", "y2"), G=2, 
                    lp=~w1+w2+w3, data=gatingsim)
```

Typically model selection issue needs to be addressed if the true model is unknown. We recommend a stepwise forwaard selection starting from a null model, i.e. fitting one bivariate gamma distribution without covariates. For example, in cases when covariates entering all gating and expert networks, a model can be fitted as

``` r
mod2 <- MBGR4(data=gatingsim, G=2, 
              l1=y1~w1+w2+w3, 
              l2=y2~w1+w2+w3, 
              l3=~w1+w2+w3, 
              l4=~w1+w2+w3, 
              lp=~w1+w2+w3)
```

Reference
---------

Hu, S., Murphy, T. B. and O'Hagan, A. (2019) Bivariate Gamma Mixture of Experts Models for Joint Claims Modeling. To appear.