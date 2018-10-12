#' Function that takes a lavaan model with standardized parameters and simulates latent scores, errors, disturbances, and observed scores
#'
#'@export
#' @param m Structural model represented by lavaan Syntax
#' @param n Number of simulated cases
#' @param ObservedOnly Return only observed data
#' @return Latent scores, errors, disturbances, and observed scores
#' @examples
#' # Lavaan model
#' m = "Latent_1 =~ 0.8 * Ob_1 + 0.8 * Ob_2"
#'
#' # simulate 10 cases
#' simStandardized(m, n = 10)
simStandardized <- function(m, n = 1000, ObservedOnly = FALSE){
  # Parameter Table
  pt <- lavParTable(m, fixed.x = F)

  # Variable Names
  vObserved <- lavNames(pt, "ov")
  vLatent <- lavNames(pt, "lv")
  vLatentExogenous <- lavNames(pt, "lv.x")
  vLatentEndogenous <- lavNames(pt, "lv.nox")
  vObservedExogenous <- lavNames(pt, "ov.x")
  vObservedEndogenous <- lavNames(pt, "ov.nox")
  if (length(vLatentEndogenous) > 0) {
    vDisturbance <- paste0("d_", vLatentEndogenous)
  } else {
    vDisturbance <- character(0)
  }

  if (length(vObservedEndogenous) > 0) {
    vError <- paste0("e_", vObservedEndogenous)
  } else {
    vError <- character(0)
  }

  # Names for A, S and new S matrices
  vA <- c(vLatentExogenous,
          vLatentEndogenous,
          vObservedExogenous,
          vObservedEndogenous)
  vS <- c(vLatentExogenous,
          vLatentEndogenous,
          vObservedExogenous,
          vObservedEndogenous)
  vnewS <- c(vLatentExogenous,
             vLatentEndogenous,
             vObservedExogenous,
             vObservedEndogenous,
             vDisturbance,
             vError)
  # Number of Variables
  k <- length(vA)

  # Initialize A matrix and exogenous correlation matrix
  ExoCor <- A <- matrix(0, k, k, dimnames = list(vA, vA))

  # Assign loadings to A
  for (i in pt[pt[, "op"] == "=~", "id"] ) {
    A[pt$rhs[i], pt$lhs[i]] <- pt$ustart[i]
  }

  # Assign regressions to A
  for (i in pt[pt[, "op"] == "~", "id"]) {
    A[pt$lhs[i], pt$rhs[i]] <- pt$ustart[i]
  }

  # Assign correlations to ExoCor
  diag(ExoCor) <- 1
  for (i in pt[pt[, "op"] == "~~", "id"]) {
    if (pt$lhs[i] != pt$rhs[i]) {
      ExoCor[pt$lhs[i], pt$rhs[i]] <- ifelse(is.na(pt$ustart[i]),
                                             0,
                                             pt$ustart[i])
      ExoCor[pt$rhs[i], pt$lhs[i]] <- ExoCor[pt$lhs[i], pt$rhs[i]]
    }

  }

  #Solving for error variances and correlation matrix

  #Column of k ones
  v1 <- matrix(1, k)

  #Initial estimate of error variances
  varS <- as.vector(v1 - (A * A) %*% v1)
  S <- diag(varS) %*% ExoCor %*% diag(varS)

  #Initial estimate of the correlation matrix
  iA <- solve(diag(k) - A)
  R <- iA %*% S %*% t(iA)

  # Set interaction count at 0
  iterations <- 0

  # Find values for S matrix
  while ((round(sum(diag(R)), 10) != k) * (iterations < 100)) {
    iA <- solve(diag(k) - A)
    R <- iA %*% S %*% t(iA)
    sdS <- diag(diag(S) ^ 0.5)
    S <- diag(diag(diag(k) - R)) + (sdS %*% ExoCor %*% sdS)
    diag(S)[diag(S) < 0] <- 0.00000001
    iterations <- iterations + 1
  }

  if (iterations  ==  100) {
    stop("Maximum iterations reached (100).")
  }

  # Assign variable names to S
  dimnames(S) <- list(vS, vS)

  # Generate data frame

  # Exogenous data
  u <- rmvnorm(n, sigma = S)
  colnames(u) <- vS

  v <- u %*% t(iA)
  #Simulated dataset
  d <- as_tibble(cbind(v,
                       u[, c(-1 * match(c(vLatentExogenous,
                                          vObservedExogenous),
                                        vS))]))
  colnames(d) <- vnewS

  dimnames(S) <- list(vS, vS)

  # Created extended A matrix
  extendA <- diag(diag(S) ^ 0.5)
  dimnames(extendA) <- list(
    c(
      vLatentExogenous,
      vLatentEndogenous,
      vObservedExogenous,
      vObservedEndogenous
    ),
    c(vLatentExogenous,
      vDisturbance,
      vObservedExogenous,
      vError)
  )
  # Remove exogenous variables
  extendA <- extendA[, c(vDisturbance, vError)]

  # bind A and extended A
  extCol <- cbind(A, extendA)
  # Append zeros so that new A will be square
  extRow <- matrix(0,
                   nrow = ncol(extendA),
                   ncol = ncol(A) + ncol(extendA))
  newA <- rbind(extCol, extRow)
  dimnames(newA) <- list(c(colnames(A), colnames(extendA)),
                         c(colnames(A), colnames(extendA)))

  # build a S matrix with 1s and 0s on the diag
  newS <- diag(c(
    rep(1, length(vLatentExogenous)),
    rep(0, length(vLatentEndogenous)),
    rep(1, length(vObservedExogenous)),
    rep(0, length(vObservedEndogenous)),
    rep(1, length(vDisturbance)),
    rep(1, length(vError))
  ))

  dimnames(newS) <- list(vnewS, vnewS)

  # Insert all off-diagonal covariances
  ExoCor <- newS[c(vLatentExogenous,
                   vObservedExogenous,
                   vDisturbance,
                   vError),
                 c(vLatentExogenous,
                   vObservedExogenous,
                   vDisturbance,
                   vError)]

  R <- solve(diag(nrow(newA)) - newA)  %*%
    newS  %*%
    t(solve(diag(nrow(newA)) - newA))

  Rxx <- R[vObserved, vObserved]
  Rxy <- R[vObserved, c(vLatent, vDisturbance, vError)]
  iRxx <- solve(Rxx)
  if (ObservedOnly) {
    d[, vObserved]
  } else {
    l <- list(
      Data = d,
      vObserved = vObserved,
      vError = vError,
      R = R,
      A = A,
      S = S,
      iterations = iterations
    )
    if (length(vLatent) > 0) {
      FScoef <- iRxx %*% Rxy
      FactorScores <- as.matrix(d[, vObserved]) %*% FScoef
      colnames(FactorScores) <- paste0(c(vLatent,
                                         vDisturbance,
                                         vError),
                                       "_FS")
      # Add factor scores to the R matrix

      #Initialise factor weight matrix
      fw <- matrix(0,
                   nrow = nrow(R),
                   ncol = ncol(FScoef),
                   dimnames = list(rownames(R), colnames(FScoef)))

      # Assign values to factor weight matrix
      fw[rownames(fw) %in% rownames(FScoef), colnames(FScoef)] <- FScoef

      # Combine factor weight matrix with identify matrix of all variables
      w <- cbind(diag(nrow(R)), fw)
      colnames(w) <- c(colnames(R), paste0(colnames(FScoef), "_FS"))

      # Correlation all variables, include factor scores
      R_all <- cov2cor(t(w) %*% R %*% w)

      # Validity coefficient (% Latent variance in factor sccores)
      FSValidity <- diag(t(FScoef) %*% Rxy)

      # Standard Errors of factor scores
      FSStandardError <- sqrt(rep(1,
                                  length(c(vLatent,
                                           vDisturbance,
                                           vError))) - FSValidity)
      paste0("se.", names(FSStandardError)) -> names(FSStandardError)

      l$Data  <- cbind(d, FactorScores)
      l$vLatent <- vLatent
      l$vDisturbance <- vDisturbance
      l$vError <- vError
      l$vFactorScores <- colnames(FactorScores)
      l$FactorScoreCoef <- FScoef
      l$R_all <- R_all
      l$R_FS <- cov2cor(t(FScoef) %*% Rxx %*% FScoef)
      l$FactorScoreValidity <- FSValidity
      l$FactorScoreSE <- FSStandardError
      l$Model <- m
    }
    l
  }
}

#' Estimate factor scores for a given profile and population model.
#'
#' @export
#' @param d observed z-scores in matrix or data.frame
#' @param m Structural model represented by lavaan Syntax
#' @return simStandardized output list with estimated factor scores from d
#' @examples
#' m = "latent_1 =~ 0.8 * ob_1 + 0.8 * ob_2 + 0.8 * ob_3"
#' d <- data.frame(ob_1 = 1, ob_2 = -0.5, ob_3 = 1.2)
#'
#' # Estimate factor scores based on this case
#' estStandardized(d = d, m = m)
estStandardized <- function(d, m){
  d_sim <- simStandardized(m = m, n = 100)
  d_FS <- as.matrix(d[, d_sim$vObserved]) %*% d_sim$FactorScoreCoef
  colnames(d_FS) <- paste0(colnames(d_FS), "_FS")
  d_sim$Data  <- cbind(d, d_FS)
  d_sim
}

#' Calculate the conditional Mahalanobis distance for any variables.
#'
#' @export
#' @param d Data.frame with the independent and dependent variables.
#' @param R Conditional correlation among variables.
#' @param Dep Vector of names of the dependent variables in your profile.
#' @param Ind Vector of names of independent variables you would like to control for.
#' @param UseFactorScores_Dep Use the factor scores for the dependent variables
#' @param UseFactorScores_Ind Use the factor scores for the independent variables
#' @param IncludeDiagnostics Return additional diagnostic information
#' @importFrom matrixcalc is.singular.matrix
#' @importFrom matrixcalc matrix.rank
#' @return conditional Mahalanobis distance, percentiles for each case based on the Chi-square distribution formed by conditional Mahalanobis distance and predicted Deps based on Inds.
#' @examples
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
#' d_demo <- simStandardized(m = m, n = 10)
#' CMahalanobis(d = d_demo$Data,
#'              R = d_demo$R_all,
#'              Dep = c("Math", "Read"),
#'              Ind = c("Gf", "Gs", "Gc"))

CMahalanobis <- function(d,
                         R,
                         Dep,
                         Ind = NULL,
                         UseFactorScores_Ind = F,
                         UseFactorScores_Dep = F,
                         IncludeDiagnostics = F) {

  if (is.list(Dep)) Dep <- unlist(Dep)
  if (is.list(Ind)) Ind <- unlist(Ind)

  # Initialize singular predictors
  Ind_Sing <- character(0)

  # Number of independent and dependent measures
  k_Ind <- length(Ind)
  k_Dep <- length(Dep)

  if (UseFactorScores_Ind) Ind_Use <- paste0(Ind, "_FS") else Ind_Use <- Ind
  if (UseFactorScores_Dep) Dep_Use <- paste0(Dep, "_FS") else Dep_Use <- Dep

  Ryy <- R[Dep, Dep]


  if (is.singular.matrix(Ryy)) {
    stop("Dependent measures are collinear. Cannot calculate the Mahalanobis Distance")
    }

  LastCondCov <- Ryy

  if (!is.null(Ind)) {
    Rxx <- R[Ind, Ind, drop = F]
    Rxy <- R[Ind, Dep, drop = F]
    Ryx <- R[Dep, Ind, drop = F]

    iRxx <- solve(Rxx)

    RegBeta <- iRxx %*% Rxy
    R2 <- colSums(RegBeta * Rxy)

    # change the name to select cases
    PredictedSubtests <- as.matrix(d[, Ind_Use]) %*% RegBeta
    SubtestDeviations <- d[, Dep_Use, drop = F] - PredictedSubtests

    #conditional variance
    CondCov <- Ryy - Ryx %*% iRxx %*% Rxy
    df <- k_Dep
    if (is.singular.matrix(CondCov)) {
      removepredictor <- function(p, Ind, Dep, R, mRank, LastCondCov) {
        if (k_Ind > 1) {
          vInd <- Ind[!(p %in% Ind)]
          Rxx <- R[vInd, vInd]
          Rxy <- R[vInd, Dep]
          Ryx <- R[Dep, vInd]
          iRxx <- solve(Rxx)
          CondCov <- Ryy - Ryx %*% iRxx %*% Rxy
          list(Remove = mRank != matrix.rank(CondCov),
               CondCov = CondCov)
        } else {
          list(Remove = TRUE,
               CondCov = LastCondCov)
        }
      }

      mRank <- matrix.rank(CondCov)
      i <- 1
      OldInd <- Ind
      while (is.singular.matrix(CondCov) | i < length(OldInd) + 1) {
          rp <- removepredictor(p  = Ind[i],
                                Ind = Ind,
                                Dep = Dep,
                                R = R,
                                mRank = mRank,
                                LastCondCov = CondCov)
          if (rp$Remove) {
            Ind <- Ind[!(Ind[i] %in% Ind)]
            CondCov <- rp$CondCov
            LastCondCov <- CondCov
            Ind_Sing <- c(Ind_Sing, Ind[i])
            df <- df - 1
          }
        i <- i + 1
      }
    }

    # Conditional Mahalanobis Distance

    dCM <- (((as.matrix(SubtestDeviations) %*%
                solve(CondCov)) * as.matrix(SubtestDeviations)) %*%
              matrix(1, nrow = k_Dep)) %>%
      sqrt %>%
      as.vector

    # Probability
    p <- pchisq(dCM ^ 2, df)

    if (IncludeDiagnostics) {
      if (k_Ind > 1) {
        # Calculate Mahalanobis distance of independent variables
        d_IndUse <- as.matrix(d[, Ind_Use, drop = F])
        dM_Ind <- (((d_IndUse %*% solve(Rxx)) * d_IndUse) %*%
          matrix(1, nrow = k_Ind)) %>%
          sqrt %>%
          as.vector
      } else dM_Ind <- NA

      # Calculate Mahalanobis distance of dependent variables
      d_Dep_Use <- as.matrix(d[, Dep_Use, drop = F])
      dM_Dep <- (((d_Dep_Use %*% solve(Ryy)) * d_Dep_Use) %*%
        matrix(1, nrow = k_Dep)) %>%
        sqrt %>%
        as.vector

      list(dCM = dCM,
           df = df,
           p = p,
           Dep = Dep_Use,
           Ind = Ind_Use,
           d_Dep = d[, Dep_Use, drop = F],
           d_Ind = d[, Ind_Use, drop = F],
           PredictedSubtests = PredictedSubtests,
           SubtestDeviations = SubtestDeviations,
           R2 = R2,
           ConditionalCovariance = CondCov,
           dM_Ind = dM_Ind,
           p_dM_Ind = pchisq(dM_Ind ^ 2, df = k_Ind),
           dM_Dep = dM_Dep,
           p_dM_Dep = pchisq(dM_Dep ^ 2, df = k_Dep))
    } else {
      list(dCM = dCM, df = df, p = p)
    }
  } else {
    d_DepUse <- as.matrix(d[, Dep_Use, drop = F])
     dCM <- (((d_DepUse %*% solve(Ryy)) * d_DepUse) %*%
               matrix(1, nrow = length(Dep))) %>%
       sqrt %>%
       as.vector
    df <- length(Dep)
    p <- pchisq(dCM ^ 2, df)
    list(dCM = dCM, df = df, p = p)
  }

}

#' Wrapper for finding out Mahalanobis distance between variables: this one gives everything for practitioners to use when they only have population relations and their clients' data
#'
#' @export
#' @param m Structural model represented by lavaan Syntax.
#' @param Dep The names of variables you would like to condition on.
#' @param Ind The names of variables of your interest.
#' @param d Profiles of interest.
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
#' maha(d = d_demo,
#'      m = m,
#'      Dep = c("Math", "Read"),
#'      Ind = c("Gc", "Gf"))
maha <- function(d, m, Dep, Ind = NULL) {
  Output <- estStandardized(d = d, m = m)
  CMahalanobis(
    Dep = Dep,
    Ind = Ind,
    R = Output$R_all,
    d = Output$Data,
    UseFactorScores_Ind = sum(Output$vLatent %in% Ind) > 1,
    UseFactorScores_Dep = sum(Output$vLatent %in% Dep) > 1
  )
}

#' Function to evaluate the accuracy of dCM with estimated factor scores
#'
#' @export
#' @param m Lavaan Syntax Object
#' @param Dep The names of variables you would like to condition on.
#' @param Ind The names of variables of your interest.
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
#' dCM_cor(m, Dep = c("Math", "Read"), Ind = c("Gc", "Gf", "Gs"), n = 100)

dCM_cor <- function(m, Dep, Ind, n = 10000) {
  tryCatch({
    # extract simulated data
    sm <- simStandardized(m, n = n)

    # get the true CMahalanobis
    dCM_true <- CMahalanobis(
      d = sm$Data,
      R = sm$R_all,
      Dep = Dep,
      Ind = Ind
    )$dCM

    # get the CMahalanobis of FS
    dCM_estimated <- CMahalanobis(
      d = sm$Data,
      R = sm$R_all,
      Dep = Dep,
      Ind = Ind,
      UseFactorScores_Ind = T,
      UseFactorScores_Dep = T
    )$dCM
    as.vector(cor(dCM_true, dCM_estimated))
  }
  , error = function(e) NA)
}

#' Confidence interval of the reliability (accuracy index)
#'
#' @export
#' @param m Population relations among variables represented by Lavaan Syntax
#' @param Dep The names of variables you would like to condition on
#' @param Ind The names of variables of your interest
#' @param replications The number of trials
#' @param sample_size The number of cases
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
#' boot_dCM(m,
#'          Dep = c("Math", "Read"),
#'          Ind = c("Gc", "Gf", "Gs"),
#'          sample_size = 100,
#'          replications = 100)
boot_dCM <- function(m,
                     Dep,
                     Ind = NULL,
                     sample_size = 10000,
                     replications = 1000) {
  out <- replicate(replications,
                   dCM_cor(
                     m = m,
                     Dep = Dep,
                     Ind = Ind,
                     n = sample_size
                   ))
  stats::quantile(out, probs = c(0.025, 0.5, 0.975))
}
