#' Calculate the conditional Mahalanobis distance for any variables.
#'
#' @export
#' @param data Data.frame with the independent and dependent variables.
#' Unless mu and sigma are specified, data are assumed to be z-scores.
#' @param R Correlation among all variables.
#' @param v_dep Vector of names of the dependent variables in your profile.
#' @param v_ind Vector of names of independent variables you would like to
#' control for.
#' @param v_ind_composites Vector of names of independent variables that are
#' composites of dependent variables
#' @param mu A vector of means. A single value means that all variables have
#' the same mean.
#' @param sigma A vector of standard deviations. A single value means that
#' all variables have the same standard deviation
#' @param use_sample_stats If TRUE, estimate R, mu, and sigma from data.
#' Only complete cases are used (i.e., no missing values in v_dep, v_ind,
#' v_ind_composites).
#' @param label optional tag for labeling output
#' @return a list with the conditional Mahalanobis distance
#' \itemize{
#' \item{\code{dCM} = Conditional Mahalanobis distance}
#' \item{\code{dCM_df} = Degrees of freedom for the conditional Mahalanobis distance}
#' \item{\code{dCM_p} = A proportion that indicates how unusual this profile is compared to profiles with the same independent variable values. For example, if \code{dCM_p} = 0.88, this profile is more unusual than 88 percent of profiles after controlling for the independent variables.}
#' \item{\code{dM_dep} = Mahalanobis distance of just the dependent variables}
#' \item{\code{dM_dep_df} = Degrees of freedom for the Mahalanobis distance of the dependent variables}
#' \item{\code{dM_dep_p} = Proportion associated with the Mahalanobis distance of the dependent variables}
#' \item{\code{dM_ind} = Mahalanobis distance of just the independent variables}
#' \item{\code{dM_ind_df} = Degrees of freedom for the Mahalanobis distance of the independent variables}
#' \item{\code{dM_ind_p} = Proportion associated with the Mahalanobis distance of the independent variables}
#' \item{\code{v_dep} = Dependent variable names}
#' \item{\code{v_ind} = Independent variable names}
#' \item{\code{v_ind_singular} = Independent variables that can be perfectly predicted from the dependent variables (e.g., composite scores)}
#' \item{\code{v_ind_nonsingular} = Independent variables that are not perfectly predicted from the dependent variables}
#' \item{\code{data} = data used in the calculations}
#' \item{\code{d_ind} = independent variable data}
#' \item{\code{d_inp_p} = Assuming normality, cumulative distribution function of the independent variables}
#' \item{\code{d_dep} = dependent variable data}
#' \item{\code{d_dep_predicted} = predicted values of the dependent variables}
#' \item{\code{d_dep_deviations = d_dep - d_dep_predicted} (i.e., residuals of the dependent variables)}
#' \item{\code{d_dep_residuals_z} = standardized residuals of the dependent variables}
#' \item{\code{d_dep_cp} = conditional proportions associated with standardized residuals}
#' \item{\code{d_dep_p} = Assuming normality, cumulative distribution function of the dependent variables}
#' \item{\code{R2} = Proportion of variance in each dependent variable explained by the independent variables}
#' \item{\code{SEE} = Standard error of the estimate for each dependent variable}
#' \item{\code{ConditionalCovariance} = Covariance matrix of the dependent variables after controlling for the independent variables}
#' \item{\code{distance_reduction = 1 - (dCM / dM_dep)} (Degree to which the independent variables decrease the Mahalanobis distance of the dependent variables. Negative reductions mean that the profile is more unusual after controlling for the independent variables. Returns 0 if \code{dM_dep} is 0.)}
#' \item{\code{variability_reduction = 1 - sum((X_dep - predicted_dep) ^ 2) / sum((X_dep - mu_dep) ^ 2)} (Degree to which the independent variables decrease the variability the dependent variables (\code{X_dep}). Negative reductions mean that the profile is more variable after controlling for the independent variables. Returns 0 if \code{X_dep == mu_dep})}
#' \item{\code{mu} = Variable means}
#' \item{\code{sigma} = Variable standard deviations}
#' \item{\code{d_person} = Data frame consisting of Mahalanobis distance data for each person}
#' \item{\code{d_variable} = Data frame consisting of variable characteristics}
#' \item{\code{label} = label slot}
#' }
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
                      mu = 0,
                      sigma = 1,
                      use_sample_stats = FALSE,
                      label = NA) {
  # Convert data to matrix
  if (is.vector(data)) {
    data <- matrix(data, nrow = 1, dimnames = list(NULL, names(data)))
  }

  data <- as.matrix(data)


  v_ind <- unique(
    c(v_ind,
      v_ind_composites))

  # Checks ----

  # Check if v_dep are in colnames of data and R
  if (!is.null(v_dep)) {
    if (!all(v_dep %in% colnames(data))) {
      v_dep_not_in_d <- paste(v_dep[!(v_ind_composites %in%
                                        colnames(data))],
                              collapse = ", ")
      stop(paste0("Some variables in v_dep are not in data: ",
                  v_dep_not_in_d))
    }

    if (!all(v_dep %in% colnames(R))) {
      v_dep_not_in_R <- paste(v_dep[!(v_dep %in%
                                        colnames(R))],
                              collapse = ", ")
      stop(paste0("Some variables in v_dep are not in R: ",
                  v_dep_not_in_R))
    }
  }



  # Check if v_ind_composites are in colnames of data
  if (!is.null(v_ind_composites)) {
    if (!all(v_ind_composites %in% colnames(data))) {
      v_ind_not_in_d <- paste(v_ind_composites[!(v_ind_composites %in%
                                                   colnames(data))],
                              collapse = ", ")
      stop(paste0("Some variables in v_ind_composites are not in data: ",
                  v_ind_not_in_d))
    }}

  # Check if v_ind_composites are in R
  if (!is.null(v_ind_composites)) {
    if (!all(v_ind_composites %in% colnames(R))) {
      v_ind_not_in_d <- paste(v_ind_composites[!(v_ind_composites %in%
                                                   colnames(R))],
                              collapse = ", ")
      stop(paste0("Some variables in v_ind_composites are not in R: ",
                  v_ind_not_in_d))
    }}

  # Check if v_ind are in colnames of data
  if (!is.null(v_ind)) {
    if (!all(v_ind %in% colnames(data))) {
      v_ind_not_in_d <- paste(v_ind[!(v_ind %in% colnames(data))],
                              collapse = ", ")
      stop(paste0("Some variables in v_ind are not in data: ",
                  v_ind_not_in_d))
    }


    if (!all(v_ind %in% colnames(R))) {
      v_ind_not_in_R <- paste(v_ind[!(v_ind %in% colnames(R))],
                              collapse = ", ")
      stop(paste0("Some variables in v_ind are not in R: ",
                  v_ind_not_in_R))
    }



    # Check if v_dep and v_ind have overlapping v
    if (any(v_dep %in% v_ind)) {
      overlap <- v_dep[v_dep %in% v_ind]
      stop(
        paste0("At least one variable is in both v_dep and v_ind: ",
               paste(overlap, collapse = " ")))
    }
  }




  # Check if v_dep are in colnames of data
  if (!all(v_dep %in% colnames(data))) {
    v_Dep_not_in_d <- paste(v_dep[!(v_dep %in% colnames(data))],
                            collapse = ", ")
    stop(paste0("Some variables in v_dep are not in data: ",
                v_Dep_not_in_d))
  }

  v_all <- c(v_ind, v_dep)
  data_k <- length(v_all)
  # Select only variables that are used
  data <- data[, v_all, drop = FALSE]

  # If data have no means specified
  if (is.null(mu) & !use_sample_stats) {
    mu <- rep(0, data_k)
    names(mu) <- colnames(data)
  } else {
    # Number of values in mu
    mu_k <- length(mu)

    # If there is only 1 value in mu, make data_k copies
    if (mu_k == 1) {
      mu <- rep(mu, data_k)
      names(mu) <- v_all
    }


    # If length of mu is unequal to the number of variables in data
    if (mu_k != data_k & mu_k > 1) {
      stop(
        paste0(
          "There are ",
          mu_k,
          " means in mu, but ",
          data_k,
          ifelse(data_k == 1, " column ", " columns "),
          "in the data. Supply 1 mean for each variable."
        )
      )
    } else {
      names(mu) <- v_all
    }
  }


  # If data have no standard deviations specified
  if (is.null(sigma) & !use_sample_stats) {
    sigma <- rep(1, data_k)
    names(sigma) <- v_all
  } else {
    # Number of values in sigma
    sigma_k <- length(sigma)

    # If there is only 1 value in sigma, make data_k copies
    if (sigma_k == 1) {
      sigma <- rep(sigma, data_k)
      names(sigma) <- v_all
    }

    # If length of sigma is unequal to the number of variables in data
    if (sigma_k != data_k & sigma_k > 1) {
      stop(
        paste0(
          "There are ",
          sigma_k,
          " means in sigma, but ",
          data_k,
          ifelse(data_k == 1, " column ", " columns "),
          "in the data. Supply 1 mean for each variable."
        )
      )
    } else {
      names(sigma) <- colnames(data)
    }
  }



  # if use_sample_stats, estimate R, mu, and sigma from data
  if (use_sample_stats) {
    if (nrow(data) < 3) stop("Sample statistics cannot be calculated with fewer than 3 rows of data. For accurate statistics, much more than 3 rows are needed.")
    d_estimate <- data[stats::complete.cases(data),v_all, drop = FALSE]
    R <- stats::cor(d_estimate)
    mu <- colMeans(d_estimate)
    sigma <- apply(d_estimate, 2, stats::sd)
  }

  # means and sd of dependent variables

  mu <- matrix(mu,
               ncol = 1,
               dimnames = list(names(mu),NULL))[v_all,, drop = FALSE]
  sigma <- matrix(sigma,
                  ncol = 1,
                  dimnames = list(names(sigma)))[v_all,, drop = FALSE]
  mu_dep <- mu[v_dep,, drop = FALSE]
  sigma_dep <- sigma[v_dep,, drop = FALSE]





  R <- R[v_all, v_all, drop = FALSE]
  # Check if R is a correlation matrix
  # Check if R has ones on the diagonal
  if (!all(diag(R) == 1)) {
    stop("R has values on its diagonal that are not ones.")}
  # Check if R is symmetric
  if (!isSymmetric(R)) stop("R is not symmetric")
  # Check if all values in R are between -1 and 1
  if (!all(R <= 1 & R >= -1)) {
    stop("Some values of R are outside the range of 1 and -1.")}


  # Make z-scores ----
  n <- nrow(data)
  ones <- matrix(rep(1, n), ncol = 1)
  colmu <- ones %*% t(mu)
  colsigma <- ones %*% t(sigma)
  d_z <- (data - colmu) / colsigma

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
  d_dep <- data[, v_dep, drop = FALSE]
  d_dep_z <- d_z[, v_dep, drop = FALSE]

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


    mu_ind <- mu[v_ind, , drop = FALSE]
    sigma_ind <- sigma[v_ind, , drop = FALSE]

    # Make matrices
    Rxx <- R[v_ind, v_ind, drop = FALSE]
    Rxy <- R[v_ind, v_dep, drop = FALSE]
    Ryx <- R[v_dep, v_ind, drop = FALSE]

    # Make sure independent covariance is nonsingluar
    if (is_singular(Rxx)) {
      stop(
        paste0(
          "Independent measures are collinear. ",
          "Cannot calculate the Mahalanobis Distance"
        )
      )
    }

    iRxx <- solve(Rxx)

    # Make regression coefficients
    reg_beta <- iRxx %*% Rxy

    # Make variance explained
    R2 <- colSums(reg_beta * Rxy)

    # Make Standard Error of Estimate
    SEE <- sqrt(1 - R2)

    # Data for just independent variables
    d_ind <- data[, v_ind, drop = FALSE]
    d_ind_z <- d_z[, v_ind, drop = FALSE]
    d_ind_p <- stats::pnorm(d_ind_z)

    # Make predicted deps
    d_predicted_z <- d_ind_z %*% reg_beta
    d_deviations_z <- d_dep_z - d_predicted_z
    d_predicted <- d_predicted_z * colsigma[, v_dep, drop = FALSE] +
      colmu[,v_dep, drop = FALSE]
    d_deviations <- d_dep - d_predicted
    variability_reduction <- 1 -  sum(d_deviations ^ 2) / sum((d_dep - mu_dep[,1, drop = TRUE]) ^ 2)


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
    dCM <- (((d_deviations_z %*%
                solve(cov_cond)) * d_deviations_z) %*%
              matrix(1, nrow = k_dep)) %>%
      sqrt %>%
      as.vector

    # Probability for Conditional Mahalanobis Distance
    dCM_p <- stats::pchisq(dCM ^ 2, dCM_df)

    # Calculate Mahalanobis distance of independent variables
    dM_ind <- (((d_ind_z %*% solve(Rxx)) * d_ind_z) %*%
                 matrix(1, nrow = k_ind)) %>%
      sqrt %>%
      as.vector

    # Probability for Mahalanobis distance of independent variables
    dM_ind_p <- stats::pchisq(dM_ind ^ 2, df = k_ind)


    # Dependent standardized residuals (z-scores)
    d_dep_residuals_z <- d_deviations_z / SEE

    # Dependent cp
    d_dep_cp <- stats::pnorm(d_dep_residuals_z)

    # Dependent p
    d_dep_p <- stats::pnorm(d_dep_z)




    d_person <- tibble::tibble(
      id = seq_len(nrow(data)),
      dCM = dCM,
      dCM_df = dCM_df,
      dCM_p = dCM_p,
      dM_dep = dM_dep,
      dM_dep_df = k_dep,
      dM_dep_p = dM_dep_p,
      dM_ind = dM_ind,
      dM_ind_df = k_ind,
      dM_ind_p = dM_ind_p
    )

    d_variable <- tibble::tibble(
      Variable = v_dep,
      mu = mu_dep[, 1],
      sigma = sigma_dep[, 1],
      R2 = R2,
      SEE = SEE) %>%
      dplyr::bind_rows(
        tibble::tibble(Variable = v_ind,
               mu = mu_ind[,1],
               sigma = sigma_ind[,1]))

    d_score <- dplyr::bind_rows(
      d_dep %>%
        tibble::as_tibble() %>%
        tibble::rowid_to_column("id") %>%
        dplyr::mutate(type = "Score"),
      d_dep_p %>%
        tibble::as_tibble() %>%
        tibble::rowid_to_column("id") %>%
        dplyr::mutate(type = "p"),
      d_dep_cp %>%
        tibble::as_tibble() %>%
        tibble::rowid_to_column("id") %>%
        dplyr::mutate(type = "cp"),
      d_predicted %>%
        tibble::as_tibble() %>%
        tibble::rowid_to_column("id") %>%
        dplyr::mutate(type = "Predicted")
      ) %>%
      dplyr::mutate(Role = "Dependent") %>%
      tidyr::pivot_longer(!!v_dep,
                          names_to = "Variable",
                          values_to = "Value") %>%
      dplyr::bind_rows(
        dplyr::bind_rows(
          d_ind %>%
            tibble::as_tibble() %>%
            tibble::rowid_to_column("id") %>%
            dplyr::mutate(type = "Score"),
          d_ind_p %>%
            tibble::as_tibble() %>%
            tibble::rowid_to_column("id") %>%
            dplyr::mutate(type = "p")) %>%
          dplyr::mutate(Role = "Independent") %>%
          tidyr::pivot_longer(!!v_ind, names_to = "Variable", values_to = "Value")
        ) %>%
      tidyr::pivot_wider(names_from = .data$type, values_from = .data$Value) %>%
      dplyr::left_join(d_variable, by = "Variable") %>%
      dplyr::mutate(Variable = factor(.data$Variable, levels = v_all))




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
      data = tibble::as_tibble(data),
      d_ind = tibble::as_tibble(d_ind),
      d_ind_p = tibble::as_tibble(d_ind_p),
      d_dep = tibble::as_tibble(d_dep),
      d_predicted = tibble::as_tibble(d_predicted),
      d_deviations = tibble::as_tibble(d_deviations),
      d_dep_residuals_z = tibble::as_tibble(d_dep_residuals_z),
      d_dep_cp = tibble::as_tibble(d_dep_cp),
      d_dep_p = tibble::as_tibble(d_dep_p),
      R2 = R2,
      SEE = SEE,
      ConditionalCovariance = cov_cond,
      distance_reduction = ifelse(dM_dep == 0, 0, 1 - (dCM / dM_dep)),
      variability_reduction = variability_reduction,
      mu = mu[, 1],
      sigma = sigma[, 1],
      d_person = d_person,
      d_score = d_score,
      d_variable = d_variable,
      label = label
    )


    class(CM) <- c("cond_maha",class(CM))
    CM



} else {
  # If there are no independent variables
  M <- list(dM_dep = dM_dep,
            dM_dep_df = k_dep,
            dM_dep_p = dM_dep_p,
            d_dep = tibble::as_tibble(d_dep),
            v_dep = v_dep,
            data = tibble::as_tibble(data),
            mu = mu[v_dep, 1],
            sigma = sigma[v_dep, 1],
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
  paste0("Conditional Mahalanobis Distance = ",
         formatC(x$dCM, 4, format = "f"),
         ", df = ",
         x$dCM_df,
         ", p = ",
         formatC(x$dCM_p, 4, format = "f"))
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
  paste0("Mahalanobis Distance = ",
         formatC(x$dM_dep, 4, format = "f"),
         ", df = ",
         x$dM_dep_df,
         ", p = ",
         formatC(x$dM_dep_p, 4, format = "f"))
}

#' Print maha class
#'
#' @param x Object of maha class
#' @noRd
#' @keywords internal
#' @export
print.maha <- function(x, ...) cat(format(x, ...))



#' Range label associated with probability
#'
#' @export
#' @param p Probability
#' @keywords internal
#' @return label string
p2label <- function(p) {
  thresholds <- purrr::map_int(p, function(x) {
    sum(x > stats::pnorm(seq(60,90,10), 100, 15)) +
      sum(x >= stats::pnorm(seq(110,140,10), 100, 15)) + 1L})
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
#' @param p probability
#' @param digits rounding digits
#'
#' @return numeric vector
#' @export
#'
#' @examples
#' proportion_round(0.01111)

proportion_round <- function(p, digits = 2) {

  lower_limit <- 0.95 * 10 ^ (-1 * digits)
  upper_limit <- 1 - lower_limit

  p1 <- dplyr::if_else(p < 0.5, p, 1 - p)

  r <- dplyr::if_else(p <= 0 | p >= 1 | p > lower_limit & p < upper_limit,
                      true = digits,
                      false = -floor(log10(abs(p1))) + (p < lower_limit))
  r <- dplyr::if_else(r > 1 & p < 1 & p > 0 & !is.na(p), r, digits)

  p2 <- stringr::str_replace(
    string = purrr::map2_chr(p, r, formatC, format = "f"),
    pattern = "^0\\.",
    replacement = ".")
  dplyr::if_else(is.na(p), NA_character_, p2)

}

#' Rounds proportions to significant digits both near 0 and 1,
#' then converts to percentiles
#'
#' @param p probability
#' @param digits rounding digits. Defaults to 2
#' @param remove_leading_zero Remove leading zero for small percentiles,
#' Defaults to TRUE
#' @param add_percent_character Append percent character. Defaults to FALSE
#' @return character vector
#' @export
#'
#' @examples
#' proportion2percentile(0.01111)

proportion2percentile <- function(p,
                                  digits = 2,
                                  remove_leading_zero = TRUE,
                                  add_percent_character = FALSE) {
  p1 <- as.character(100 * as.numeric(proportion_round(as.numeric(p),
                                                       digits = digits)))
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
#' @param x The results of the cond_maha function.
#' @param ... Arguments passed to print function
#' @param p_tail The proportion of the tail to shade
#' @param family Font family.
#' @param score_digits Number of digits to round scores.
#' @importFrom rlang .data
#'
plot.cond_maha <- function(x,
                           ...,
                           p_tail = 0,
                           family = "serif",
                           score_digits = ifelse(min(x$sigma) >= 10, 0, 2)) {

  if (length(unique(x$d_score$id)) > 1) stop("Can only plot one case at a time")

  break_width <- max(x$sigma)
  break_min <- min(x$mu - 10 * x$sigma)
  break_max <- max(x$mu + 10 * x$sigma)
  minor_break_width <- ifelse(break_width %% 3 == 0, break_width / 3, break_width / 2)
  major_breaks <- seq(break_min, break_max, break_width)
  minor_breaks <- seq(break_min, break_max, minor_break_width)

  label_independent <- paste0("list(Independent~italic(d[M])==\"",
                              formatC(x$dM_ind,
                                   digits = 2,
                                   format = "f"),
                              "\",italic(p)==\"",
                              proportion_round(x$dM_ind_p),
                              "\")")


  label_dependent <- paste0("list(Dependent~italic(d[M])==\"",
                            formatC(x$dM_dep,
                                    digits = 2,
                                    format = "f"),
                            "\",italic(p)==\"",
                            proportion_round(x$dM_dep_p),
                            "\")")

  x$d_score %>%
    dplyr::mutate(SD = ifelse(test = is.na(.data$SEE),
                              yes = .data$sigma,
                              no = .data$SEE * .data$sigma),
           yhat = ifelse(is.na(.data$Predicted), .data$mu, .data$Predicted),
           id = factor(.data$id),
           Role = factor(.data$Role,
                         levels = c("Independent", "Dependent"),
                         labels = c(label_independent, label_dependent))) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$Variable,
                                 y = .data$Score,
                                 fill = .data$Role)) +
    ggplot2::facet_grid(cols = ggplot2::vars(!!quote(Role)),
               scales = "free",
               space = "free",
               labeller = ggplot2::label_parsed) +
    ggnormalviolin::geom_normalviolin(
      mapping = ggplot2::aes(
        mu = .data$yhat,
        sigma = .data$SD,
        face_right = .data$Role == label_dependent,
        face_left = .data$Role != label_dependent),
      p_tail = p_tail,
      fill = "gray68",
      width = 0.85) +
    ggnormalviolin::geom_normalviolin(
      mapping = ggplot2::aes(
        mu = .data$mu,
        sigma = .data$sigma,
        face_right = .data$Role != label_dependent),
      fill = "gray88",
      width = 0.85,
      p_tail = p_tail) +
    ggplot2::geom_point(mapping = ggplot2::aes(color = .data$id)) +
    ggplot2::geom_text(
      mapping = ggplot2::aes(
        label = formatC(.data$Score, score_digits, format = "f"),
        color = .data$id
      ),
      vjust = -0.5,
      family = family
    ) +
    ggplot2::geom_text(
      mapping = ggplot2::aes(
        color = .data$id,
        label = dplyr::if_else(
          .data$Role == label_independent,
          "",
          paste0(
            "italic(c*p)=='",
            proportion_round(.data$cp),
            "'"))),
      vjust = 2.3,
      size = 3,
      parse = TRUE,
      family = family) +
    ggplot2::geom_text(
      mapping = ggplot2::aes(
        color = .data$id,
        label = paste0(
          "italic(phantom(c)*p)=='",
            proportion_round(.data$p),
          "'"
        )
      ),
      vjust = 1.3,
      size = 3,
      parse = TRUE,
      family = family) +
    ggplot2::scale_y_continuous("Scores",
                                breaks = major_breaks,
                                minor_breaks = minor_breaks) +
    ggplot2::scale_x_discrete(NULL,
                     expand = ggplot2::expansion(add = 1)) +
    ggplot2::labs(title = bquote(list(
      Conditional ~ Mahalanobis ~ Distance ~ (italic(d[CM])) == .(formatC(x$dCM, 2, format = "f")),
      italic(p) == .(proportion_round(x$dCM_p))
    )),
    subtitle = bquote(list(
      Mahalanobis ~ Distance ~ Reduction == .(paste0(round(100 * x$distance_reduction), "%")),
      Euclidean ~ Distance ~ Reduction == .(paste0(round(100 * x$variability_reduction), "%"))
    )),
    caption = expression(
      list(
        italic(p) == "Population proportion",
        italic(c * p) == "Conditional proportion",
        italic(d[M]) == "Mahalanobis Distance"
      )
    )) +
    ggplot2::theme_light(base_family = family) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_color_grey()
}

#' Plot objects of the maha class (i.e, the results of the cond_maha function
#' using dependent variables only).
#'
#' @export
#' @param x The results of the cond_maha function.
#' @param ... Arguments passed to print function
#' @param p_tail Proportion in violin tail (defaults to 0).
#' @param family Font family.
#' @param score_digits Number of digits to round scores.
#' @importFrom rlang .data
#'
plot.maha <- function(x,
                      ...,
                      p_tail = 0,
                      family = "serif",
                      score_digits = ifelse(min(x$sigma) >= 10, 0, 2)) {

  if (nrow(x$d_dep) > 1) stop("Can only plot one case at a time")

  break_width <- max(x$sigma)
  break_min <- min(x$mu - 10 * x$sigma)
  break_max <- max(x$mu + 10 * x$sigma)
  minor_break_width <- ifelse(break_width %% 3 == 0, break_width / 3, break_width / 2)
  major_breaks <- seq(break_min, break_max, break_width)
  minor_breaks <- seq(break_min, break_max, minor_break_width)

  d_stats <- tibble::tibble(Variable = x$v_dep, mu = x$mu, sigma = x$sigma)

  x$d_dep %>%
    tibble::rownames_to_column("id") %>%
    dplyr::mutate(id = factor(.data$id)) %>%
    tidyr::pivot_longer(-.data$id,
                        names_to = "Variable",
                        values_to = "Score") %>%
    dplyr::left_join(d_stats, by = "Variable") %>%
    dplyr::mutate(
      z = (.data$Score - .data$mu) / .data$sigma,
      z_p = stats::pnorm(.data$z),
      in_tail = .data$z_p < p_tail / 2 | .data$z_p > (1 - p_tail / 2)) %>%
    ggplot2::ggplot(ggplot2::aes(.data$Variable, .data$Score)) +
    ggnormalviolin::geom_normalviolin(data = d_stats,
                                      ggplot2::aes(x = .data$Variable,
                                          mu = .data$mu,
                                          sigma = .data$sigma),
                                      inherit.aes = FALSE,
                                      fill = "gray65",
                                      p_tail = p_tail) +
    ggplot2::geom_point(mapping = ggplot2::aes(shape = .data$in_tail)) +
    ggplot2::geom_text(
      mapping = ggplot2::aes(
        label = formatC(.data$Score, score_digits, format = "f"),
        color = .data$id
      ),
      vjust = -0.5,
      family = family
    ) +
    ggplot2::geom_text(
      mapping = ggplot2::aes(
        label = paste0("italic(p)=='",
                       proportion_round(.data$z_p),
                       "'")),
      vjust = 1.3,
      size = 3,
      parse = TRUE,
      family = family) +
    ggplot2::theme_minimal(base_family = family) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_color_grey() +
    ggplot2::scale_shape_manual(values = c(16, 0)) +
    ggplot2::scale_y_continuous("Scores",
                                breaks = major_breaks,
                                minor_breaks = minor_breaks) +
    ggplot2::scale_x_discrete(NULL,
                              expand = ggplot2::expansion(add = 1)) +
    ggplot2::labs(title = bquote(list(
      Mahalanobis == .(formatC(x$dM_dep, 2, format = "f")),
      italic(p) == .(proportion_round(x$dM_dep_p))
    )),
    caption = expression(
      list(
        italic(p) == "Population proportion"
      )
    ))
  }
