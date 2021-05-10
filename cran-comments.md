## Test environments

-   local OS X install, R 3.5.1
-   ubuntu 12.04 (on travis-ci), R 3.5.1
-   win-builder (devel and release)

## R CMD check results

Duration: 36.4s

0 errors v \| 0 warnings v \| 0 notes v

R CMD check succeeded

## Addressing concerns from a previous build

From Uwe Ligges:

> "The unusualprofile packages calculates" is not grammatical and redundant, hence please omit it. Rather explain the methodology. Is there some reference about the method you can add in the Description field in the form Authors (year) <doi:.....>?


Response:

The description field now reads:

> Calculates a Mahalanobis distance for every row of a set of outcome variables (Mahalanobis, 1936 <doi: 10.1007/s13171-019-00164-5>). The conditional Mahalanobis distance is calculated using a conditional covariance matrix (i.e., a covariance matrix of the outcome variables after controlling for a set of predictors). Plotting the output of the cond_maha function can help identify which elements of a profile are unusual after controlling for the predictors.
