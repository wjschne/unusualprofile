library(testthat)
library(simstandard)
# lavaan model with three indicators of a latent variable
model <- "
X =~ 0.7 * X_1 + 0.5 * X_2 + 0.8 * X_3
Y =~ 0.7 * Y_1 + 0.5 * Y_2 + 0.8 * Y_3
Y ~ 0.6 * X
"

# Randomly generated case
d <- sim_standardized(
  model,
  n = 1,
  observed = TRUE,
  latent = FALSE,
  errors = FALSE,
  composites = TRUE) %>%
  rename(X = "X_Composite",
         Y = "Y_Composite")

# Model-implied correlation matrix
R <- sim_standardized_matrices(model)$Correlations$R_all
v_dep = c("X_1", "X_2", "X_3",
          "Y_1", "Y_2", "Y_3")
v_ind_composites = c("X", "Y")

test_that("Example works", {
  expect_silent(cond_maha(d, R,v_dep = v_dep, v_ind_composites = v_ind_composites))
})

test_that("R is a correlation matrix", {
  R1 <- R
  R1[1,1] <- 1.1
  expect_error(cond_maha(d, R1,v_dep = v_dep, v_ind_composites = v_ind_composites))
  R1[1,1] <- -1.1
  expect_error(cond_maha(d, R1,v_dep = v_dep, v_ind_composites = v_ind_composites))
  R1[1,1] <- 0.5
  expect_error(cond_maha(d, R1,v_dep = v_dep, v_ind_composites = v_ind_composites))
  R1[1,1] <- 1
  R1[1,2] <- 0.5
  R1[2,1] <- 0.6
  expect_error(cond_maha(d, R1,v_dep = v_dep, v_ind_composites = v_ind_composites))
})

test_that("v_ind, v_ind_composites, and v_dep are in colnames of data", {
  expect_error(cond_maha(d, R,v_dep = v_dep, v_ind <- "fred", v_ind_composites = v_ind_composites))
  expect_error(cond_maha(d, R,v_dep = c(v_dep,"fred"), v_ind_composites = v_ind_composites))
  expect_error(cond_maha(d, R,v_dep = v_dep, v_ind_composites = c(v_ind_composites, "fred")))
})

test_that("v_ind and v_dep have no overlapping names", {
  expect_error(cond_maha(d, R,v_dep = v_dep, v_ind <- "X_1", v_ind_composites = v_ind_composites))
  expect_error(cond_maha(d, R,v_dep = c(v_dep,"X"), v_ind_composites = v_ind_composites))
  expect_error(cond_maha(d, R,v_dep = v_dep, v_ind_composites = c(v_ind_composites, "X_1")))
})

