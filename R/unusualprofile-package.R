#' unusualprofile: Calculates Conditional Mahalanobis Distances
#'
#' @docType package
#' @name unusualprofile
#'
#' @description
#'
#' The unusualprofile package calculates the unusualness of score profiles conditioned on
#' a set of predictor variables
#'
#'
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom mvtnorm rmvnorm
#' @importFrom tibble as_tibble
#' @importFrom tibble rowid_to_column
#' @importFrom stats cov2cor
#' @importFrom stats pchisq
#' @importFrom stats rbeta
#' @importFrom stats cor
#' @importFrom tidyr unite
#' @importFrom glue glue
#' @importFrom matrixcalc is.singular.matrix
#' @importFrom matrixcalc matrix.rank
#' @import lavaan
#' @import dplyr
#' @import purrr
#' @author Feng Ji
#' @author W. Joel Schneider
#'
NULL
