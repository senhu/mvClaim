
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

Reference
---------

Hu, S., Murphy, T. B. and O'Hagan, A. (2019) Bivariate Gamma Mixture of Experts Models for Joint Claims Modeling. To appear.
