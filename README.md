
<!-- README.md is generated from README.Rmd. Please edit that file -->

# unusualprofile

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/ggnormalviolin)](https://cran.r-project.org/package=unusualprofile)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Travis build
status](https://travis-ci.org/wjschne/unusualprofile.svg?branch=master)](https://travis-ci.org/wjschne/unusualprofile)
[![Codecov test
coverage](https://codecov.io/gh/wjschne/unusualprofile/branch/master/graph/badge.svg)](https://codecov.io/gh/wjschne/unusualprofile?branch=master)
<!-- badges: end -->

The goal of unusualprofile is to calculate conditional Mahalanobis
distances and related statistics. Such statistics can help find cases
with unusual profiles of multivariate normal data.

# Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("wjschne/unusualprofile")
```

# Example

To use the unusualprofile package, all that is needed is to know the
correlations, means, and standard deviations among a set of continuous
variables and at least one row of data from that set of variables.

## Preliminary work: Generate a model and case data

This section has nothing to with the unusualprofile package per se. For
now, we need to create some data to work with. Suppose we have set of
variables that have the following relationships:

![Multivariate normal model](man/figures/rm_model.svg)

First, we load the packages that we will use later:

``` r
library(unusualprofile)
library(simstandard)
library(ggnormalviolin)
library(dplyr)
```

Using the [simstandard
package](https://wjschne.github.io/simstandard/articles/simstandard_tutorial.html),
we generate 1 case from a standardized multivariate normal data set. We
use [lavaan syntax](http://lavaan.ugent.be/tutorial/syntax1.html) to
match the structure in the figure:

``` r
# lavaan model with three indicators of a latent variable
model <- "
X =~ 0.7 * X_1 + 0.5 * X_2 + 0.8 * X_3
Y =~ 0.8 * Y_1 + 0.7 * Y_2 + 0.9 * Y_3
Y ~ 0.6 * X
"
```

Now we generate a case from this multivariate distribution:

``` r
# Ensure that random data will be the same as in the example.
set.seed(281)

# Randomly generated case
d <- sim_standardized(
  model, 
  n = 1, 
  observed = TRUE, 
  latent = FALSE, 
  errors = FALSE, 
  composites = TRUE) 
```

The case in the `d` variable looks like this:

| X\_1 | X\_2 | X\_3 | Y\_1 |   Y\_2 | Y\_3 | X\_Composite | Y\_Composite |
| ---: | ---: | ---: | ---: | -----: | ---: | -----------: | -----------: |
|    2 | 1.45 |  2.7 | 1.66 | \-0.84 | 1.81 |         2.59 |         1.01 |

The model-implied correlation matrix:

``` r

# Model-implied correlation matrix
R <- sim_standardized_matrices(model)$Correlations$R_all
```

|              | X\_1 | X\_2 | X\_3 | Y\_1 | Y\_2 | Y\_3 | X\_Composite | Y\_Composite |
| ------------ | ---: | ---: | ---: | ---: | ---: | ---: | -----------: | -----------: |
| X\_1         | 1.00 | 0.35 | 0.56 | 0.34 | 0.29 | 0.38 |         0.81 |         0.39 |
| X\_2         | 0.35 | 1.00 | 0.40 | 0.24 | 0.21 | 0.27 |         0.74 |         0.28 |
| X\_3         | 0.56 | 0.40 | 1.00 | 0.38 | 0.34 | 0.43 |         0.83 |         0.44 |
| Y\_1         | 0.34 | 0.24 | 0.38 | 1.00 | 0.56 | 0.72 |         0.40 |         0.87 |
| Y\_2         | 0.29 | 0.21 | 0.34 | 0.56 | 1.00 | 0.63 |         0.35 |         0.84 |
| Y\_3         | 0.38 | 0.27 | 0.43 | 0.72 | 0.63 | 1.00 |         0.46 |         0.90 |
| X\_Composite | 0.81 | 0.74 | 0.83 | 0.40 | 0.35 | 0.46 |         1.00 |         0.47 |
| Y\_Composite | 0.39 | 0.28 | 0.44 | 0.87 | 0.84 | 0.90 |         0.47 |         1.00 |

## Using the `cond_maha` function

Now we specify the correlations (`R`), means (`mu`), standard deviations
(`sigma`). independent variables (`v_ind`), and dependent variables
(`v_dep`). In this case, the indpendent variables are composite scores
summarizing the dependent variables.

``` r
# Conditional Mahalanobis distance
cm <- cond_maha(data = d, 
          R = R,
          mu = 0,
          sigma = 1,
          v_ind_composites = c("X_Composite", "Y_Composite"),
          v_dep = c("X_1", "X_2", "X_3",
                    "Y_1", "Y_2", "Y_3"))

cm
#> Conditional Mahalanobis Distance = 3.1117, df = 4, p = 0.9539

# Plot
plot_cond_maha(cm)
```

<img src="man/figures/README-example-1.svg" width="100%" />
