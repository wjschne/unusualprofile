#' Wrapper for finding out Mahalanobis distance between variables:
#' this one gives everything for practitioners to use when they only have
#' population relations and their clients' data
#'
#' @export
#' @param data Profiles of interest.
#' @param model Structural model represented by lavaan Syntax.
#' @param v_dep The names of variables you would like to condition on.
#' @param v_ind The names of variables of your interest.
#' @importFrom rlang .data
#' @return conditional Mahalanobis distance, percentiles for each case
#' based on the Chi-square distribution formed by conditional Mahalanobis
#' distance and predicted Deps based on Inds.
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

  # Select only observed scores,
  # Add factor scores,
  # Rename them to latent variables
  d <- data[, sm$v_names$v_observed, drop = FALSE] %>%
    simstandard::add_factor_scores(model) %>%
    dplyr::rename_with(stringr::str_remove, pattern = "_FS")

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

  # Calculate Conditional Mahalanobis distances
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
    dplyr::mutate(
      v_dep = purrr::map(.data$value, "v_dep"),
      v_ind = purrr::map(.data$value, "v_ind"),
      label = purrr::map(.data$value, "label")
    ) %>%
    dplyr::select(-.data$value, -.data$Type) %>%
    purrr::pmap(.f = cond_maha,
                d = data,
                R = sm$Correlations$R_all)


  tibble::tibble(
    Measure = purrr::map_chr(cm_all, "label"),
    dCM = purrr::map_dbl(cm_all, "dCM"),
    `dCM Percentile` = purrr::map_dbl(cm_all, "dCM_p"),
    `dCM Level` = p2label(.data$`dCM Percentile`),
    dM = purrr::map_dbl(cm_all, "dM_dep"),
    `dM Percentile` = purrr::map_dbl(cm_all, "dM_dep_p"),
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
#' @return Correlation between the conditional Mahalanobis distance
#' calculated by using the true scores and the conditional Mahalanobis
#' calculated by using estimated factor scores
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

    R_all <- simstandard::sim_standardized_matrices(model)$Correlations$R_all

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
      v_dep = paste0(v_dep, "_FS"),
      v_ind = paste0(v_ind, "_FS")
    )$dCM
    stats::cor(dCM_true, dCM_estimated)

  }
  , error = function(e) NA)
}

#' Confidence interval of the reliability (accuracy index)
#'
#' @export
#' @param model character string representing a model in lavaan syntax
#' @param v_dep Vector of names of the dependent variables in your profile.
#' @param v_ind Vector of names of independent variables you would like
#' to control for.
#' @param v_ind_composites Vector of names of independent variables
#' that are composites of dependent variables
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
                     v_ind_composites = NULL,
                     sample_size = 1000,
                     replications = 1000) {
  if (is.null(v_ind)) {
    v_ind_FS <- NULL
  } else {
    v_ind_FS <- paste0(v_ind, "_FS")
  }

  if (is.null(v_ind_composites)) {
    v_ind_composites_FS <- NULL
  } else {
    v_ind_composites_FS <- paste0(
      v_ind_composites,
      "_FS")
  }

  if (is.null(v_dep)) {
    v_dep_FS <- NULL
  } else {
    v_dep_FS <- paste0(v_dep, "_FS")
  }

  sm <- simstandard::sim_standardized(
    m = model,
    n = sample_size * replications,
    factor_scores = TRUE,
    errors = FALSE, observed = FALSE)

  R_all <- simstandard::sim_standardized_matrices(model)$Correlations$R_all

  tidyr::crossing(id = 1:sample_size, sample_id = 1:replications) %>%
    dplyr::bind_cols(sm) %>%
    dplyr::group_by(.data$sample_id) %>%
    tidyr::nest() %>%
    dplyr::mutate(dcm_cor = purrr::map_dbl(.data$data, function(d) {
      dcm_obs <- cond_maha(
        data = d %>%
          dplyr::select(!!unique(c(v_ind_FS,
                                   v_dep_FS))),
        R = R_all,
        v_dep = v_dep_FS,
        v_ind = v_ind_FS,
      )$dCM

      dcm_lat <- cond_maha(
        data = d %>%
          dplyr::select(!!unique(c(v_ind,
                                   v_dep))),
        R = R_all,
        v_dep = v_dep,
        v_ind = v_ind
      )$dCM

      stats::cor(dcm_obs, dcm_lat)

    })
    ) %>%
    dplyr::pull(dcm_cor) %>%
    stats::quantile(probs = c(0.025, 0.5, 0.975))



}
