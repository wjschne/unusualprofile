#functions for validation

#' Automatically calculate the accuracy of one condition.
#' @export
#' @param Dep The names of variables you would like to condition on.
#' @param Ind  The names of variables of your interest.
#' @param d  Dataset created by using simStandarized.
#' @return Correlation between the conditional Mahalanobis distance calculated by using the true scores and the conditional Mahalanobis calculated by using estimated factor scores
#' @examples
#' SimModel <- "
#' Gc =~ 0.85 * Gc1 + 0.68 * Gc2 + 0.8 * Gc3
#' Gf =~ 0.8 * Gf1 + 0.9 * Gf2 + 0.8 * Gf3
#' Gs =~ 0.7 * Gs1 + 0.8 * Gs2 + 0.8 * Gs3
#' Read =~ 0.66 * Read1 + 0.85 * Read2 + 0.91 * Read3
#' Math =~ 0.4 * Math1 + 0.9 * Math2 + 0.7 * Math3
#' Gc ~ 0.6 * Gf + 0.1 * Gs
#' Gf ~ 0.5 * Gs
#' Read ~ 0.4 * Gc + 0.1 * Gf
#' Math ~ 0.2 * Gc + 0.3 * Gf + 0.1 * Gs"
#' d_demo <- simStandardized(SimModel, 10)
#' CM_cor(Dep = c("Math", "Read"),Ind = c("Gc", "Gf", "Gs") , d_demo)

CM_cor <- function(Dep, Ind = NULL, d){

  # extract simulated data
  R <- d[["R"]]
  cor_data <- d[["Data"]]

  # get the true CMahalanobis
  TrueCM <- CMahalanobis(Dep = Dep, Ind = Ind, R = R, d = cor_data)

  # get the CMahalanobis of FS
  EstCM <- CMahalanobis_FS(Dep = Dep, Ind = Ind, R = R, d = cor_data)
  TrueCM <- TrueCM[["dCM"]]
  EstCM <- EstCM[["dCM"]]

  # calculate the reliability
  cor <- cor(TrueCM, EstCM)
  cor
}

#' Function that removes all fixed values in a lavaan model
#' @export
#' @param m Population relations among variables represented by Lavaan Syntax.
#' @return Lavaan object without specified coefficients
#' @examples
#' SimModel <- "
#' Gc =~ 0.85 * Gc1 + 0.68 * Gc2 + 0.8 * Gc3
#' Gf =~ 0.8 * Gf1 + 0.9 * Gf2 + 0.8 * Gf3
#' Gs =~ 0.7 * Gs1 + 0.8 * Gs2 + 0.8 * Gs3
#' Read =~ 0.66 * Read1 + 0.85 * Read2 + 0.91 * Read3
#' Math =~ 0.4 * Math1 + 0.9 * Math2 + 0.7 * Math3
#' Gc ~ 0.6 * Gf + 0.1 * Gs
#' Gf ~ 0.5 * Gs
#' Read ~ 0.4 * Gc + 0.1 * Gf
#' Math ~ 0.2 * Gc + 0.3 * Gf + 0.1 * Gs"
#' Sim2Free(SimModel)

Sim2Free <- function(m){
  m %>%
    lavaanify(fixed.x = F) %>%
    filter(.data$lhs != .data$rhs) %>%
    group_by(.data$lhs, .data$op) %>%
    summarise(rhs = paste(.data$rhs, collapse = " + ")) %>%
    arrange(desc(.data$op)) %>%
    tidyr::unite("l", .data$lhs, .data$op, .data$rhs, sep = " ") %>%
    pull(.data$l) %>%
    paste(collapse = "\n")
}

#' Function to make a lavaan formula into a model for simulation that takes a vector of parameters as input
#' @export
#' @param SimModel Population relations among variables represented by Lavaan Syntax.
#' @return Lavaan object without specified coefficients
#' @examples
#' SimModel <- "
#' Gc =~ 0.85 * Gc1 + 0.68 * Gc2 + 0.8 * Gc3
#' Gf =~ 0.8 * Gf1 + 0.9 * Gf2 + 0.8 * Gf3
#' Gs =~ 0.7 * Gs1 + 0.8 * Gs2 + 0.8 * Gs3
#' Read =~ 0.66 * Read1 + 0.85 * Read2 + 0.91 * Read3
#' Math =~ 0.4 * Math1 + 0.9 * Math2 + 0.7 * Math3
#' Gc ~ 0.6 * Gf + 0.1 * Gs
#' Gf ~ 0.5 * Gs
#' Read ~ 0.4 * Gc + 0.1 * Gf
#' Math ~ 0.2 * Gc + 0.3 * Gf + 0.1 * Gs"
#' sim2glue(SimModel)

sim2glue <- function(SimModel){
  SimModel %>%
    lavaanify(fixed.x = F) %>%
    filter(.data$lhs != .data$rhs) %>%
    rowid_to_column(var = "ID") %>%
    mutate(ustart = paste0("{para",
                           "[",
                           sprintf( "%02d", .data$ID ),
                           "]",
                           "} * ")) %>%
    group_by(.data$op, .data$lhs) %>%
    summarise(rhs = paste(.data$ustart, .data$rhs, collapse = " + ")) %>%
    arrange(desc(.data$op)) %>%
    arrange(.data$rhs) %>%
    unite("l", .data$lhs, .data$op, .data$rhs, sep = " ") %>%
    pull(.data$l) %>%
    paste(collapse = "\n")
}

#' Function to stuck parameters to the lavaan model
#' @export
#' @param SimModel Population relations among variables represented by Lavaan Syntax.
#' @param para Parameters for simulation studies.
#' @return A list of lavaan object with specified parameters.

glue2simf <- function(SimModel, para){
  glue(
    sim2glue(SimModel)
  )
}


#' Create a tibble to hold everything for simulation
#' @export
#' @param SimModel Population relations among variables represented by Lavaan Syntax
#' @param dPar Parameters we would like to test in simulation studies
#' @param n The number of cases
#' @return all parameters, a lavaan object with specified parameters and data generated for each condition by calling simStandarized
#' @examples
#' SimModel <- "
#' Gc =~ 0.85 * Gc1 + 0.68 * Gc2 + 0.8 * Gc3
#' Gf =~ 0.8 * Gf1 + 0.9 * Gf2 + 0.8 * Gf3
#' Gs =~ 0.7 * Gs1 + 0.8 * Gs2 + 0.8 * Gs3
#' Read =~ 0.66 * Read1 + 0.85 * Read2 + 0.91 * Read3
#' Math =~ 0.4 * Math1 + 0.9 * Math2 + 0.7 * Math3
#' Gc ~ 0.6 * Gf + 0.1 * Gs
#' Gf ~ 0.5 * Gs
#' Read ~ 0.4 * Gc + 0.1 * Gf
#' Math ~ 0.2 * Gc + 0.3 * Gf + 0.1 * Gs"
#' beta_measure1 <- makePara(1000, 9000, n = 10, k = 15)
#' beta_structure1 <- makePara(1000, 1000, n = 10, k = 8)
#' cbind(beta_measure1, beta_structure1)-> par_beta1
#' colnames(par_beta1) <- paste0("x", 1:23)
#' simMaha(par_beta1, SimModel, n = 10)

simMaha <- function(dPar, SimModel, n = 100){
  by_row(dPar, glue2simf, SimModel = SimModel) %>%
    mutate(data = map(.data$.out,
                      simStandardized,
                      n = n))
}

#' Simulate data and cor between true conditional Mahalanobis distance and the estimated
#' @export
#' @param SimModel Population relations among variables represented by Lavaan Syntax
#' @param dPar Parameters we would like to test in simulation studies
#' @param Dep The names of variables you would like to condition on
#' @param Ind  The names of variables of your interest
#' @param n The number of cases
#' @return all parameters, lavaan object with specified parameters and data generated for each condition by calling simStandarized and the accuracy index for each condition
#' @examples
#' SimModel <- "
#' Gc =~ 0.85 * Gc1 + 0.68 * Gc2 + 0.8 * Gc3
#' Gf =~ 0.8 * Gf1 + 0.9 * Gf2 + 0.8 * Gf3
#' Gs =~ 0.7 * Gs1 + 0.8 * Gs2 + 0.8 * Gs3
#' Read =~ 0.66 * Read1 + 0.85 * Read2 + 0.91 * Read3
#' Math =~ 0.4 * Math1 + 0.9 * Math2 + 0.7 * Math3
#' Gc ~ 0.6 * Gf + 0.1 * Gs
#' Gf ~ 0.5 * Gs
#' Read ~ 0.4 * Gc + 0.1 * Gf
#' Math ~ 0.2 * Gc + 0.3 * Gf + 0.1 * Gs"
#' beta_measure1 <- makePara(1000, 9000, n = 10, k = 15)
#' beta_structure1 <- makePara(1000, 1000, n = 10, k = 8)
#' par_beta1 <- cbind(beta_measure1, beta_structure1)
#' colnames(par_beta1) <- paste0("x", 1:23)
#' simCor(dPar = par_beta1, SimModel, Dep = c("Read", "Math"),Ind = c("Gc", "Gs", "Gf"), n = 10)

simCor <- function(dPar, SimModel, Dep, Ind = NULL, n = 100){
  dPar %>%
    simMaha(SimModel, n) %>%
    mutate(cor = map_dbl(.data$data,
                         CM_cor,
                         Dep = Dep,
                         Ind = Ind))
}

#' function to make a set of parameters that follows beta distribution
#' @export
#' @param from starting point for the parameters
#' @param to ending point for the parameters
#' @param n number of cases
#' @param k  number of variables.
#' @return all parameters, lavaan object with specified parameters and data generated for each condition by calling simStandarized and the accuracy index for each condition
#' @examples
#' makePara(1000, 9000, n = 10, k = 15)

makePara <- function(from, to, n, k){
  list1 <- map(seq(from = from,
                   to = to,
                   length.out = n),
               function(p) rbeta(k, p, 10000 - p))
  data.frame(t(sapply(list1, c)))
}

#' confidence interval of the reliability (accuracy index)
#' @export
#' @param SimModel Population relations among variables represented by Lavaan Syntax
#' @param Dep The names of variables you would like to condition on
#' @param Ind  The names of variables of your interest
#' @param size  The number of trials
#' @param n  The number of cases
#' @return simulated 95% confidence interval
#' @examples
#' SimModel <- "
#' Gc =~ 0.85 * Gc1 + 0.68 * Gc2 + 0.8 * Gc3
#' Gf =~ 0.8 * Gf1 + 0.9 * Gf2 + 0.8 * Gf3
#' Gs =~ 0.7 * Gs1 + 0.8 * Gs2 + 0.8 * Gs3
#' Read =~ 0.66 * Read1 + 0.85 * Read2 + 0.91 * Read3
#' Math =~ 0.4 * Math1 + 0.9 * Math2 + 0.7 * Math3
#' Gc ~ 0.6 * Gf + 0.1 * Gs
#' Gf ~ 0.5 * Gs
#' Read ~ 0.4 * Gc + 0.1 * Gf
#' Math ~ 0.2 * Gc + 0.3 * Gf + 0.1 * Gs"
#' boot(SimModel, Dep = c("Math", "Read"),Ind = c("Gc", "Gf", "Gs"), size = 100, n = 100)

boot <- function(SimModel, Dep, Ind = NULL, size = 100, n = 100){
  con_cor <- function(SimModel, Dep, Ind, size){
    d <- simStandardized(SimModel, size)
    CM_cor(Dep, Ind, d)
  }
  out <- replicate(n, con_cor(SimModel, Dep, Ind, size))
  stats::quantile(out, probs = c(0.025, 0.975))
}
