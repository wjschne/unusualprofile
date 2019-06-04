
<!-- README.md is generated from README.Rmd. Please edit that file -->

# unusualprofile

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/ggnormalviolin)](https://cran.r-project.org/package=unusualprofile)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Travis build
status](https://travis-ci.org/wjschne/unusualprofile.svg?branch=master)](https://travis-ci.org/wjschne/unusualprofile)
<!-- badges: end -->

The goal of unusualprofile is to calculate conditional Mahalanobis
distances and related statistics. Such statistics can help find cases
with unusual profiles of multivariate normal data.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("wjschne/unusualprofile")
```

## Example

Using the simstandard package, we generate 1 case from a standardized
multivariate normal data set.

``` r
library(unusualprofile)
library(simstandard)
library(ggnormalviolin)
library(tidyverse)
library(extrafont)
set.seed(48)

# lavaan model with three indicators of a latent variable
model <- "
X =~ 0.7 * X_1 + 0.5 * X_2 + 0.8 * X_3
Y =~ 0.7 * Y_1 + 0.5 * Y_2 + 0.8 * Y_3
Y ~ 0.6 * X
"

# Randomly generated case
d <- sim_standardized(
  model, 
  n = 1, 
  observed = TRUE, 
  latent = FALSE, 
  errors = FALSE, 
  composites = TRUE) %>% 
  rename(X = "X_Composite",
         Y = "Y_Composite")

# Model-implied correlation matrix
R <- sim_standardized_matrices(model)$Correlations$R_all
```

Now we specify the correlations, independent variables, and dependent
variables. In this case, the indpendent variables are composite scores
summarizing the dependent variables.

``` r
# Conditional Mahalanobis distance
cm <- cond_maha(data = d, 
          R = R,
          v_dep = c("X_1", "X_2", "X_3", 
                    "Y_1", "Y_2", "Y_3"),
          v_ind_composites = c("X", "Y"))

cm
#> Conditional Mahalanobis Distance = 3.0834, df = 4, p = 0.9504

# Plot
plot_cond_maha(cm)
```

<img src="man/figures/README-unnamed-chunk-2-1.svg" width="100%" />
