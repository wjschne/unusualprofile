#' Calculate the conditional Mahalanobis distance for any variables.
#'
#' @export
#' @param data Data.frame with the independent and dependent variables.
#' @param R Correlation among all variables.
#' @param v_dep Vector of names of the dependent variables in your profile.
#' @param v_ind Vector of names of independent variables you would like to control for.
#' @param label optional tag for labeling output
#' @return conditional Mahalanobis distance, percentiles for each case based on the Chi-square distribution formed by conditional Mahalanobis distance and predicted Deps based on Inds.
#' @examples
#' library(unusualprofile)
#' library(simstandard)
#'
#' m <- "
#' Gc =~ 0.85 * Gc1 + 0.68 * Gc2 + 0.8 * Gc3
#' Gf =~ 0.8 * Gf1 + 0.9 * Gf2 + 0.8 * Gf3
#' Gs =~ 0.7 * Gs1 + 0.8 * Gs2 + 0.8 * Gs3
#' Read =~ 0.66 * Read1 + 0.85 * Read2 + 0.91 * Read3
#' Math =~ 0.4 * Math1 + 0.9 * Math2 + 0.7 * Math3
#' Gc ~ 0.6 * Gf + 0.1 * Gs
#' Gf ~ 0.5 * Gs
#' Read ~ 0.4 * Gc + 0.1 * Gf
#' Math ~ 0.2 * Gc + 0.3 * Gf + 0.1 * Gs"
#' # Generate 10 cases
#' d_demo <- simstandard::sim_standardized(m = m, n = 10)
#'
#' # Get model-implied correlation matrix
#' R_all <- simstandard::sim_standardized_matrices(m)$Correlations$R_all
#'
#' cond_maha(data = d_demo,
#'           R = R_all,
#'           v_dep = c("Math", "Read"),
#'           v_ind = c("Gf", "Gs", "Gc"))

cond_maha <- function(data,
                      R,
                      v_dep,
                      v_ind = NULL,
                      label = NA) {
  # Convert d to matrix
  if (is.vector(data))
    data <- matrix(data, nrow = 1)
  data <- as.matrix(data)

  # Check if v_ind are in colnames d
  if (!all(v_ind %in% colnames(data))) {
    v_Ind_not_in_d <-
      paste(v_ind[!(v_ind %in% colnames(data))], collapse = ", ")
    stop(paste0("Some variables in v_ind are not in d: ", v_Ind_not_in_d))
  }

  # Check if v_dep are in colnames d
  if (!all(v_dep %in% colnames(data))) {
    v_Dep_not_in_d <-
      paste(v_dep[!(v_dep %in% colnames(data))], collapse = ", ")
    stop(paste0("Some variables in v_dep are not in d: ", v_Dep_not_in_d))
  }

  # Select covariance among dependent variables
  Ryy <- R[v_dep, v_dep, drop = FALSE]

  # Make sure dependent covariance is nonsingluar
  if (is_singular(Ryy)) {
    stop(
      paste0(
        "Dependent measures are collinear. ",
        "Cannot calculate the Mahalanobis Distance"
      )
    )
  }


  # Data for just dependent variables
  d_dep <- data[, v_dep, drop = F]

  # Number of dependent variables
  k_dep <- length(v_dep)

  # (Unconditional) Mahalanobis distance of dependent variables
  dM_dep <- (((d_dep %*% solve(Ryy)) * d_dep) %*%
               matrix(1, nrow = k_dep)) %>%
    sqrt %>%
    as.vector
  # Probability for Mahalanobis distance of dependent variables
  dM_dep_p <- stats::pchisq(dM_dep ^ 2, df = k_dep)

  if (!is.null(v_ind)) {
    # If there are independent variables,
    # calculate conditional Mahalanobis distance.

    # Make matrices
    Rxx <- R[v_ind, v_ind, drop = FALSE]
    Rxy <- R[v_ind, v_dep, drop = FALSE]
    Ryx <- R[v_dep, v_ind, drop = FALSE]
    iRxx <- solve(Rxx)

    # Make regression coefficients
    reg_beta <- iRxx %*% Rxy

    # Make variance explained
    R2 <- colSums(reg_beta * Rxy)

    # Make Standard Error of Estimate
    SEE <- sqrt(1 - R2)

    # Data for just independent variables
    d_ind <- data[, v_ind, drop = F]

    # Make predicted deps
    d_predicted <- d_ind %*% reg_beta
    d_deviations <- d_dep - d_predicted

    # Conditional Variance
    cov_cond <- Ryy - Ryx %*% iRxx %*% Rxy
    k_ind <- length(v_ind)
    dCM_df <- k_dep

    # Initialize predictor vectors
    v_ind_singular <- character(0)
    v_ind_nonsingular <- character(0)

    if (is_singular(cov_cond)) {
      # Function to see if adding a predictor makes the
      # conditional covariance singular
      addpredictor <- function(p, oldInd, v_dep, R) {
        vInd <- c(p, oldInd)
        Rxx <- R[vInd, vInd]
        Rxy <- R[vInd, v_dep]
        Ryx <- R[v_dep, vInd]
        iRxx <- solve(Rxx)
        cov_cond <- Ryy - Ryx %*% iRxx %*% Rxy
        !is_singular(cov_cond)
      }

      # Add predictors that do not make conditional
      # covariance singular
      for (i in v_ind) {
        if (addpredictor(i, v_ind_nonsingular, v_dep, R)) {
          v_ind_nonsingular <- c(v_ind_nonsingular, i)
        } else {
          v_ind_singular <- c(v_ind_singular, i)
          dCM_df <- dCM_df - 1
        }
      }

      # Make final conditional covariance
      if (length(v_ind_nonsingular) > 0) {
        Rxx <- R[v_ind_nonsingular, v_ind_nonsingular]
        Rxy <- R[v_ind_nonsingular, v_dep]
        Ryx <- R[v_dep, v_ind_nonsingular]
        iRxx <- solve(Rxx)
        cov_cond <- Ryy - Ryx %*% iRxx %*% Rxy
      } else
        cov_cond <- Ryy
    }


    # Conditional Mahalanobis Distance
    dCM <- (((d_deviations %*%
                solve(cov_cond)) * d_deviations) %*%
              matrix(1, nrow = k_dep)) %>%
      sqrt %>%
      as.vector

    # Probability for Conditional Mahalanobis Distance
    dCM_p <- stats::pchisq(dCM ^ 2, dCM_df)

    # Calculate Mahalanobis distance of independent variables
    dM_ind <- (((d_ind %*% solve(Rxx)) * d_ind) %*%
                 matrix(1, nrow = k_ind)) %>%
      sqrt %>%
      as.vector

    # Probability for Mahalanobis distance of independent variables
    dM_ind_p <- stats::pchisq(dM_ind ^ 2, df = k_ind)


    # Dependent standardized residuals (z-scores)
    d_dep_residuals_z <- d_deviations / SEE

    # Dependent p
    d_dep_residuals_p <- stats::pnorm(d_dep_residuals_z)


    CM <- list(
      dCM = dCM,
      dCM_df = dCM_df,
      dCM_p = dCM_p,
      dM_dep = dM_dep,
      dM_dep_df = k_dep,
      dM_dep_p = dM_dep_p,
      dM_ind = dM_ind,
      dM_ind_df = k_ind,
      dM_ind_p = dM_ind_p,
      v_dep = v_dep,
      v_ind = v_ind,
      v_ind_singular = v_ind_singular,
      v_ind_nonsingular = v_ind_nonsingular,
      d_dep = tibble::as_tibble(d_dep),
      d_ind = tibble::as_tibble(d_ind),
      d_predicted = tibble::as_tibble(d_predicted),
      d_deviations = tibble::as_tibble(d_deviations),
      d_dep_residuals_z = tibble::as_tibble(d_dep_residuals_z),
      d_dep_residuals_p = tibble::as_tibble(d_dep_residuals_p),
      R2 = R2,
      SEE = SEE,
      ConditionalCovariance = cov_cond,
      variability_explained = 1 - (dCM / dM_dep) ^ 2,
      label = label
    )

    class(CM) <- c("cond_maha",class(CM))
    CM



} else {
  # If there are no independent variables
  M <- list(dM_dep = dM_dep,
       dM_dep_df = k_dep,
       dM_dep_p = dM_dep_p,
       label = label)
  class(M) <- c("maha",class(M))
  M
  }
}

#' Format cond_maha class
#'
#' @param x Object of cond_maha class
#' @noRd
#' @keywords internal
#' @export
format.cond_maha <- function(x, ...) {
  paste0("Conditional Mahalanobis Distance = ",formatC(x$dCM, 4, format = "f"), ", df = ", x$dCM_df, ", p = ", formatC(x$dCM_p, 4, format = "f"))
}

#' Print cond_maha class
#'
#' @param x Object of cond_maha class
#' @noRd
#' @keywords internal
#' @export
print.cond_maha <- function(x, ...) cat(format(x, ...), "\n")


#' Format maha class
#'
#' @param x Object of maha class
#' @noRd
#' @keywords internal
#' @export
format.maha <- function(x, ...) {
  paste0("Mahalanobis Distance = ",formatC(x$dM_dep, 4, format = "f"), ", df = ", x$dM_dep_df, ", p = ", formatC(x$dM_dep_p, 4, format = "f"))
}

#' Print maha class
#'
#' @param x Object of maha class
#' @noRd
#' @keywords internal
#' @export
print.maha <- function(x, ...) cat(format(x, ...), "\n")



#' Wrapper for finding out Mahalanobis distance between variables: this one gives everything for practitioners to use when they only have population relations and their clients' data
#'
#' @export
#' @param data Profiles of interest.
#' @param model Structural model represented by lavaan Syntax.
#' @param v_dep The names of variables you would like to condition on.
#' @param v_ind The names of variables of your interest.
#' @importFrom rlang .data
#' @return conditional Mahalanobis distance, percentiles for each case based on the Chi-square distribution formed by conditional Mahalanobis distance and predicted Deps based on Inds.
#' @examples
#' # Standardized structural model in lavaan syntax
#' m <- "
#' Gc =~ 0.85 * Gc_1 + 0.68 * Gc_2 + 0.80 * Gc_3
#' Gf =~ 0.80 * Gf_1 + 0.90 * Gf_2 + 0.80 * Gf_3
#' Read =~ 0.66 * Read_1 + 0.85 * Read_2 + 0.91 * Read_3
#' Math =~ 0.40 * Math_1 + 0.90 * Math_2 + 0.70 * Math_3
#' Gc ~ 0.60 * Gf
#' Read ~ 0.40 * Gc + 0.10 * Gf
#' Math ~ 0.20 * Gc + 0.30 * Gf
#' "
#'
#' # Put observed scores in data.frame
#' d_demo <- data.frame(
#'           Gc_1 = -1,
#'           Gc_2 = 0.5,
#'           Gc_3 = -0.2,
#'           Gf_1 = 1.1,
#'           Gf_2 = 1.3,
#'           Gf_3 = 2,
#'           Read_1 = -0.5,
#'           Read_2 = -1,
#'           Read_3 = -1.4,
#'           Math_1 = 1.1,
#'           Math_2 = 1.3,
#'           Math_3 = 0.7
#' )
#' unusualness(
#'      data = d_demo,
#'      model = m,
#'      v_dep = c("Math", "Read"),
#'      v_ind = c("Gc", "Gf"))
unusualness <- function(
  data,
  model,
  v_dep,
  v_ind
  ) {

  if (nrow(data) > 1) stop("This function works only with 1 row of data")

  # Model matrices
  sm <- simstandard::sim_standardized_matrices(model)

  data <- simstandard::add_factor_scores(data, model)

  # Observed variables
  v_dep_obs <- lavaan::lavaanify(model) %>%
    dplyr::filter(.data$op == "=~") %>%
    dplyr::filter(.data$lhs %in% v_dep) %>%
    dplyr::pull(.data$rhs) %>%
    unique()

  v_ind_obs <- lavaan::lavaanify(model) %>%
    dplyr::filter(.data$op == "=~") %>%
    dplyr::filter(.data$lhs %in% v_ind) %>%
    dplyr::pull(.data$rhs) %>%
    unique()

  # Calculate Conditional mahalanobis distances
  cm_all <- list(
    ind_lat2dep_lat = list(
      v_dep = v_dep,
      v_ind = v_ind,
      label = "Latent Dependent $\\leftarrow$ Latent Independent"),
    ind_obs2dep_obs = list(
      v_dep = v_dep_obs,
      v_ind = v_ind_obs,
      label = "Observed Dependent $\\leftarrow$ Observed Independent"),
    dep_lat2dep_obs = list(
      v_dep = v_dep_obs,
      v_ind = v_dep,
      label = "Observed Dependent $\\leftarrow$ Latent Dependent"),
    ind_lat2ind_obs = list(
      v_dep = v_ind_obs,
      v_ind = v_ind,
      label = "Observed Independent $\\leftarrow$ Latent Independent"),
    lat2obs = list(
      v_dep = c(v_ind_obs, v_dep_obs),
      v_ind = c(v_ind, v_dep),
      label = "All Observed $\\leftarrow$ All Latent")
    ) %>%
    tibble::enframe(name = "Type") %>%
    dplyr::mutate(v_dep = purrr::map(.data$value, "v_dep"),
           v_ind = purrr::map(.data$value, "v_ind"),
           label = purrr::map(.data$value, "label")) %>%
    dplyr::select(-.data$value,-.data$Type) %>%
    purrr::pmap(.f = cond_maha,
                d = data,
                R = sm$Correlations$R_all)


  tibble::tibble(
    Measure = purrr::map_chr(cm_all,"label"),
    dCM = purrr::map_dbl(cm_all,"dCM"),
    `dCM Percentile` = purrr::map_dbl(cm_all,"dCM_p"),
    `dCM Level` = p2label(.data$`dCM Percentile`),
    dM = purrr::map_dbl(cm_all,"dM_dep"),
    `dM Percentile` = purrr::map_dbl(cm_all,"dM_dep_p"),
    `dM Level` = p2label(.data$`dM Percentile`)
  )




}

#' Function to evaluate the accuracy of dCM with estimated factor scores
#'
#' @export
#' @param model Lavaan Syntax Object
#' @param v_dep The names of variables you would like to condition on.
#' @param v_ind The names of variables of your interest.
#' @param n Sample size of simulated data
#' @return Correlation between the conditional Mahalanobis distance calculated by using the true scores and the conditional Mahalanobis calculated by using estimated factor scores
#' @examples
#' m <- "
#' Gc =~ 0.85 * Gc_1 + 0.68 * Gc_2 + 0.80 * Gc_3
#' Gf =~ 0.80 * Gf_1 + 0.90 * Gf_2 + 0.80 * Gf_3
#' Gs =~ 0.70 * Gs_1 + 0.80 * Gs_2 + 0.80 * Gs_3
#' Read =~ 0.66 * Read_1 + 0.85 * Read_2 + 0.91 * Read_3
#' Math =~ 0.40 * Math_1 + 0.90 * Math_2 + 0.70 * Math_3
#' Gc ~ 0.60 * Gf + 0.10 * Gs
#' Gf ~ 0.50 * Gs
#' Read ~ 0.40 * Gc + 0.10 * Gf
#' Math ~ 0.20 * Gc + 0.30 * Gf + 0.10 * Gs
#' "
#' dcm_cor(m, v_dep = c("Math", "Read"), v_ind = c("Gc", "Gf", "Gs"), n = 100)

dcm_cor <- function(model, v_dep, v_ind, n = 10000) {
  tryCatch({
    # extract simulated data
    sm <- simstandard::sim_standardized(
      m = model,
      n = n,
      factor_scores = TRUE)

    R_all  = simstandard::sim_standardized_matrices(model)$Correlations$R_all

    # get the true CMahalanobis
    dCM_true <- cond_maha(
      data = sm,
      R = R_all,
      v_dep = v_dep,
      v_ind = v_ind
    )$dCM

    # get the CMahalanobis of FS
    dCM_estimated <- cond_maha(
      data = sm,
      R = R_all,
      v_dep = paste0(v_dep,"_FS"),
      v_ind = paste0(v_ind,"_FS")
    )$dCM
    stats::cor(dCM_true, dCM_estimated)

  }
  , error = function(e) NA)
}

#' Confidence interval of the reliability (accuracy index)
#'
#' @export
#' @param model Population relations among variables represented by Lavaan Syntax
#' @param v_dep The names of variables you would like to condition on
#' @param v_ind The names of variables of your interest
#' @param replications The number of trials
#' @param sample_size The number of cases
#' @importFrom rlang .data
#' @return simulated 95% confidence interval
#' @examples
#' m <- "
#' Gc =~ 0.85 * Gc_1 + 0.68 * Gc_2 + 0.80 * Gc_3
#' Gf =~ 0.80 * Gf_1 + 0.90 * Gf_2 + 0.80 * Gf_3
#' Gs =~ 0.70 * Gs_1 + 0.80 * Gs_2 + 0.80 * Gs_3
#' Read =~ 0.66 * Read_1 + 0.85 * Read_2 + 0.91 * Read_3
#' Math =~ 0.40 * Math_1 + 0.90 * Math_2 + 0.70 * Math3
#' Gc ~ 0.60 * Gf + 0.10 * Gs
#' Gf ~ 0.50 * Gs
#' Read ~ 0.40 * Gc + 0.10 * Gf
#' Math ~ 0.20 * Gc + 0.30 * Gf + 0.10 * Gs
#' "
#' boot_dcm(model = m,
#'          v_dep = c("Math", "Read"),
#'          v_ind = c("Gc", "Gf", "Gs"),
#'          sample_size = 100,
#'          replications = 100)
boot_dcm <- function(model,
                     v_dep,
                     v_ind = NULL,
                     sample_size = 1000,
                     replications = 1000) {


  sm <- simstandard::sim_standardized(
    m = model,
    n = sample_size * replications,
    factor_scores = TRUE)

  R_all  = simstandard::sim_standardized_matrices(model)$Correlations$R_all

  tidyr::crossing(id = 1:sample_size, sample_id = 1:replications) %>%
    dplyr::bind_cols(sm) %>%
    dplyr::group_by(.data$sample_id) %>%
    tidyr::nest() %>%
    dplyr::mutate(dcm_cor = purrr::map_dbl(.data$data, function(d) {
      dcm_obs <- cond_maha(
      data = d,
      R = R_all,
      v_dep = v_dep,
      v_ind = v_ind
      )$dCM

      dcm_lat <- cond_maha(
        data = d,
        R = R_all,
        v_dep = paste0(v_dep,"_FS"),
        v_ind = paste0(v_ind,"_FS")
      )$dCM

      stats::cor(dcm_obs, dcm_lat)

      })
      ) %>%
    dplyr::pull(dcm_cor) %>%
    stats::quantile(probs = c(0.025, 0.5, 0.975))



}

#' Range label associated with probability
#'
#' @export
#' @param p Probability
#' @keywords internal
#' @return label string
p2label <- function(p) {
  thresholds <- purrr::map_int(p, function(x) sum(x > stats::pnorm(seq(60,140,10)[-5], 100, 15)) + 1L)
  vlabels <- c(
    "Extremely Low",
    "Very Low",
    "Low",
    "Low Average",
    "Average",
    "High Average",
    "High",
    "Very High",
    "Extremely High")
  vlabels[thresholds]
}


#' Test if matrix is singular
#'
#' @export
#' @param x matrix
#' @keywords internal
#' @return logical
is_singular <- function(x) {
  det(x) < .Machine$double.eps
}
