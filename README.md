---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# clustTMB

<!-- badges: start -->
[![R-CMD-check](https://github.com/Andrea-Havron/clustTMB/workflows/R-CMD-check/badge.svg)](https://github.com/Andrea-Havron/clustTMB/actions)
[![codecov](https://codecov.io/gh/Andrea-Havron/clustTMB/branch/main/graph/badge.svg?token=kfHYcSdj88)](https://codecov.io/gh/Andrea-Havron/clustTMB)
<!-- badges: end -->

clustTMB fits finite mixture models that incorporate both covariate and random effect
structure within either the gating, the expert, or both components of the model. The
clustTMB model is built in Template Model Builder (TMB) and implemented in R using formula
syntax similar to other finite mixture modeling packages such as mclust. The TMB structure
of the finite mixture model allows for fast computations of random effects by using the
Laplace Approximation to calculate the marginal likelihood of fixed effects after integrating
out random components. Spatial random effects are estimated through approximations of
the Gaussian Markov Random Fields via the FEM-SPDE method introduced by INLA.

## Installation

You can install the development version from [GitHub](https://github.com/) with:


```r
remotes::install_github("Andrea-Havron/clustTMB", build_vignettes = TRUE)
```

The clustTMB package relies on functions from R-INLA, which is not available on CRAN. Details on installing R-INLA are available on [this](https://www.r-inla.org/download-install) website. To install the stable version, use:

```r
install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
```


## Example

A basic example is illustrated with the faithful dataset:


```r
library(clustTMB)
#> Error in library(clustTMB): there is no package called 'clustTMB'
data("faithful")
m1 <- clustTMB(response = faithful, covariance.structure = "VVV")
#> Error in clustTMB(response = faithful, covariance.structure = "VVV"): could not find function "clustTMB"
```



```
#> Error in clustTMB(response = faithful, covariance.structure = "VVV"): could not find function "clustTMB"
#> Error in eval(expr, envir, enclos): object 'm1' not found
```

<div class="figure">
<img src="man/figures/README-pressure-1.png" alt="plot of chunk pressure" width="100%" />
<p class="caption">plot of chunk pressure</p>
</div>

For more detailed examples, see vignettes:

```r
## Not Run

# install.packages("utils")
utils::browseVignettes("clustTMB")

## End Not Run
```
