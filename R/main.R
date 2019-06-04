#' Calculate the conditional Mahalanobis distance for any variables.
#'
#' @export
#' @param data Data.frame with the independent and dependent variables. Unless mu and sigma are specified, data are assumed to be z-scores.
#' @param R Correlation among all variables.
#' @param v_dep Vector of names of the dependent variables in your profile.
#' @param v_ind Vector of names of independent variables you would like to control for.
#' @param v_ind_composites Vector of names of independent variables that are composites of dependent variables
#' @param mu A vector of means. A single value means that all variables have the same mean.
#' @param sigma A vector of standard deviations. A single value means that all variables have the same standard deviation
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
                      v_ind_composites = NULL,
                      mu = NULL,
                      sigma = NULL,
                      label = NA) {
  # Convert d to matrix
  if (is.vector(data))
    data <- matrix(data, nrow = 1)
  data <- as.matrix(data)


  v_ind <- unique(c(v_ind, v_ind_composites))

  # Check if v_ind are in colnames of data
  if (!is.null(v_ind)) {
    if (!all(v_ind %in% colnames(data))) {
      v_ind_not_in_d <-
        paste(v_ind[!(v_ind %in% colnames(data))], collapse = ", ")
      stop(paste0("Some variables in v_ind are not in d: ", v_ind_not_in_d))
    }

    # Check if v_dep and v_ind have overlapping v
    if (any(v_dep %in% v_ind)) {
      overlap <- v_dep[v_dep %in% v_ind]
      stop(paste0("At least one variable is in both v_dep and v_ind: ", paste(overlap, collapse = " ")))
    }
  }




  # Check if v_dep are in colnames d
  if (!all(v_dep %in% colnames(data))) {
    v_Dep_not_in_d <-
      paste(v_dep[!(v_dep %in% colnames(data))], collapse = ", ")
    stop(paste0("Some variables in v_dep are not in d: ", v_Dep_not_in_d))
  }

  v_all <- c(v_ind, v_dep)

  data_k <- ncol(data)

  # If data have no means specified
  if (is.null(mu)) {
    mu <- rep(0, data_k)
  } else {
    # Number of values in mu
    mu_k <- length(mu)

    # If there is only 1 value in mu, make data_k copies
    if (mu_k == 1)
      mu <- rep(mu, data_k)

    # If length of mu is unequal to the number of variables in data
    if (mu_k != data_k & mu_k > 1) {
      stop(
        paste0(
          "There are ",
          mu_k,
          " means in mu, but ",
          data_n,
          ifelse(data_k == 1, " column ", " columns "),
          "in the data. Supply 1 mean for each variable."
        )
      )
    }
  }


  # If data have no standard deviations specified
  if (is.null(sigma)) {
    sigma <- rep(0, data_k)
  } else {
    # Number of values in sigma
    sigma_k <- length(sigma)

    # If there is only 1 value in sigma, make data_k copies
    if (sigma_k == 1)
      sigma <- rep(sigma, data_k)

    # If length of sigma is unequal to the number of variables in data
    if (sigma_k != data_k & sigma_k > 1) {
      stop(
        paste0(
          "There are ",
          sigma_k,
          " means in sigma, but ",
          data_n,
          ifelse(data_k == 1, " column ", " columns "),
          "in the data. Supply 1 mean for each variable."
        )
      )
    }
  }

  # Check if R is a correlation matrix
  # Check if R has ones on the diagonal
  if (!all(diag(R) == 1)) stop("R has values on its diagonal that are not ones.")
  # Check if R is symmetric
  if (!isSymmetric(R)) stop("R is not symmetric")
  # Check if all values in R are between -1 and 1
  if (!all(R <= 1 & R >= -1)) stop("Some values of R are outside the range of 1 and -1.")


  # Make z-scores
  n <- nrow(data)
  ones <- matrix(rep(1, n), ncol = 1)
  colmu <- ones %*% matrix(mu, nrow = 1)
  colsd <- ones %*% matrix(sigma, nrow = 1)
  d_z <- (data - colmu) / colsd

  d_z <- d_z[, v_all, drop = FALSE]
  data <- data[, v_all, drop = FALSE]

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
  d_dep_z <- d_z[, v_dep, drop = F]

  # Number of dependent variables
  k_dep <- length(v_dep)

  # (Unconditional) Mahalanobis distance of dependent variables
  dM_dep <- (((d_dep_z %*% solve(Ryy)) * d_dep_z) %*%
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
    if (is.null(v_ind_composites)) v_ind_composites <- character(0)
    v_ind_singular <- v_ind_composites
    v_ind_nonsingular <- character(0)
    v_ind_try <- setdiff(v_ind, v_ind_singular)
    dCM_df <- dCM_df - length(v_ind_singular)

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
      for (i in v_ind_try) {
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
print.maha <- function(x, ...) format(x, ...)



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

#' Rounds proportions to significant digits both near 0 and 1
#'
#' @param p probabiity
#' @param digits rounding digits
#'
#' @return numeric vector
#' @export
#'
#' @examples
#' proportion_round(0.01111)

proportion_round <- function(p, digits = 2) {
  p1 <- round(p, digits)
  lower_limit <- 0.95 * 10 ^ (-1 * digits)
  upper_limit <- 1 - lower_limit
  p1[p > upper_limit & p <= 1] <- 1 - signif(1 - p[p > upper_limit & p <= 1], digits - 1)
  p1[p < lower_limit & p >= 0] <- signif(p[p < lower_limit & p >= 0], digits - 1)
  p1
}

#' Rounds proportions to significant digits both near 0 and 1, then converts to percentiles
#'
#' @param p probabiity
#' @param digits rounding digits
#'
#' @return chracter vector
#' @export
#'
#' @examples
#' prop2percentile(0.01111)

proportion2percentile <- function(p,
                                  digits = 2,
                                  remove_leading_zero = TRUE,
                                  add_percent_character = FALSE) {
  p1 <- as.character(100 * proportion_round(p, digits = digits))
  if (remove_leading_zero) {
    p1 <- gsub(pattern = "^0\\.",
               replacement = ".",
               x = p1)
  }


  if (add_percent_character) {
    p1 <- paste0(p1,"%")
  }

  p1

}


#' Plot the variables from the results of the cond_maha function.
#'
#' @export
#' @param cm The results of the cond_maha function.
#' @param family Font family.
plot_cond_maha <- function(cm, family = "serif") {
  d <- bind_rows(
    cm$d_ind %>%
      rowid_to_column("id") %>%
      gather("Variable", "Score",-id) %>%
      mutate(
        Predicted = 0,
        SD = 1,
        p = pnorm(Score, 0, 1),
        Role = "Independent"
      ),
    cm$d_dep %>%
      rowid_to_column("id") %>%
      gather("Variable", "Score",-id) %>%
      left_join(
        cm$d_predicted %>%
          rowid_to_column("id") %>%
          gather("Variable", "Predicted",-id),
        by = c("Variable", "id")
      ) %>%
      left_join(tibble(Variable = names(cm$SEE),
                       SD = cm$SEE),
                by = "Variable") %>%
      left_join(
        cm$d_dep_residuals_p %>%
          rowid_to_column("id") %>%
          gather("Variable", "p",-id),
        by = c("Variable", "id")
      ) %>%
      mutate(Role = "Dependent")) %>%
    mutate(Role = fct_inorder(Role),
           id = factor(id),
           p = proportion_round(p))

    ggplot(d, aes(Variable, Score, fill = Role)) +
        facet_grid(
          cols = vars(Role),
          scales = "free",
          space = "free"
        ) +
        ggnormalviolin::geom_normalviolin(
          aes(
            mu = Predicted,
            sigma = SD,
            face_right = Role == "Dependent",
            face_left = Role != "Dependent"
          ),
          fill = "gray90"
        ) +
        ggnormalviolin::geom_normalviolin(aes(
          mu = 0,
          sigma = 1,
          face_right = Role != "Dependent"
        ),
        fill = "gray65") +
        geom_point(mapping = aes(color = id)) +
    geom_text(
      mapping = aes(label = formatC(Score, 2, format = "f"),
                    color = id),
      vjust = -0.5,
      family = family
    ) +
    geom_text(
      mapping = aes(
        color = id,
        label = if_else(
          Role == "Independent",
          "",
          paste0(
            "italic(c*p)=='",
            str_replace_all(
              p,
              pattern = "0\\.",
              replacement = "."
            ),
            "'"
          )
        )
      ),
      vjust = 2.3,
      size = 3,
      parse = TRUE,
      family = family
    ) +
    geom_text(
      mapping = aes(
        color = id,
        label = paste0(
          "italic(p)=='",
          str_replace_all(
            proportion_round(pnorm(Score, 0, 1)),
            "0\\.",
            "."),
          "'")
      ),
      vjust = 1.3,
      size = 3,
      parse = TRUE,
      family = family
    ) +
    scale_y_continuous("z-Scores") +
    scale_x_discrete(NULL,
                     expand = expand_scale(add = 1)) +
    labs(title = bquote(list(
      Conditional ~ Mahalanobis == .(formatC(cm$dCM, 2, format = "f")),
      italic(p) == .(proportion_round(cm$dCM_p))
    )),
    caption = expression(
      list(
        italic(p) == "Population proportion",
        italic(c * p) == "Conditional proportion"
      )
    )) +
    theme_light(base_family = family) +
    theme(legend.position = "none")
}

