#' An example data.frame
#'
#' A dataset with 1 row of data for a single case.
#'
#' @format A data frame with 1 row and 8 variables:
#' \describe{
#' \item{X_1}{A predictor variable}
#' \item{X_2}{A predictor variable}
#' \item{X_3}{A predictor variable}
#' \item{Y_1}{An outcome variable}
#' \item{Y_2}{An outcome variable}
#' \item{Y_3}{An outcome variable}
#' \item{X}{A latent predictor variable}
#' \item{Y}{A latent outcome variable}
#' }
"d_example"

#' An example correlation matrix
#'
#' A correlation matrix used for demonstration purposes
#' It is the model-implied correlation matrix for this structural model:
#' `X =~ 0.7 * X_1 + 0.5 * X_2 + 0.8 * X_3`
#' `Y =~ 0.8 * Y_1 + 0.7 * Y_2 + 0.9 * Y_3`
#' `Y ~ 0.6 * X`
#'
#' @format A matrix with 8 rows and 8 columns:
#' \describe{
#' \item{X_1}{A predictor variable}
#' \item{X_2}{A predictor variable}
#' \item{X_3}{A predictor variable}
#' \item{Y_1}{An outcome variable}
#' \item{Y_2}{An outcome variable}
#' \item{Y_3}{An outcome variable}
#' \item{X}{A latent predictor variable}
#' \item{Y}{A latent outcome variable}
#' }
"R_example"
