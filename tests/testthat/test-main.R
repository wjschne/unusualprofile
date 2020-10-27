# Make data ----

library(testthat)
library(simstandard)
library(dplyr)
library(unusualprofile)
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
  latent = TRUE,
  errors = FALSE,
  composites = TRUE
)

# Model-implied correlation matrix
v_dep <- c("X_1", "X_2", "X_3",
          "Y_1", "Y_2", "Y_3")
v_ind_composites <- c("X_Composite", "Y_Composite")
v_ind <- NULL
v_all <- c(v_ind, v_ind_composites, v_dep)
R <- sim_standardized_matrices(model)$Correlations$R_all[v_all,
                                                         v_all,
                                                         drop = FALSE]

# Works  and Input ----
test_that("Example works", {
  expect_silent(cond_maha(
    data = d,
    R = R,
    v_dep = v_dep,
    v_ind_composites = v_ind_composites
  ))
})

test_that("data as matrix okay", {
  expect_silent(cond_maha(
    data = as.matrix(d),
    R = R,
    v_dep = v_dep,
    v_ind_composites = v_ind_composites
  ))
})

test_that("data as vector okay", {
  expect_silent(cond_maha(
    as.matrix(d)[1, ],
    R = R,
    v_dep = v_dep,
    v_ind_composites = v_ind_composites
  ))
})

test_that("R is a correlation matrix", {
  R1 <- R
  R1[1, 2] <- 1.1
  R1[2, 1] <- 1.1
  expect_error(
    cond_maha(
      data = d,
      R = R1,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    ),
    "Some values of R are outside the range of 1 and -1."
  )
  R1 <- R
  R1[1, 2] <- -1.1
  R1[2, 1] <- -1.1
  expect_error(
    cond_maha(
      data = d,
      R = R1,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    ),
    "Some values of R are outside the range of 1 and -1."
  )
  R1 <- R
  R1[1, 1] <- 0.5
  expect_error(
    cond_maha(
      data = d,
      R = R1,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    ),
    "R has values on its diagonal that are not ones."
  )
  R1 <- R
  R1[1, 2] <- 0.5
  R1[2, 1] <- 0.6
  expect_error(
    cond_maha(
      data = d,
      R = R1,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    ),
    "R is not symmetric"
  )
})

# Variable names ----

test_that("v_ind, v_ind_composites, and v_dep are in colnames of data", {
  expect_error(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind = "fred",
      v_ind_composites = v_ind_composites
    ),
    "Some variables in v_ind are not in data"
  )
  expect_error(
    cond_maha(
      data = d,
      R = R,
      v_dep = c(v_dep, "fred"),
      v_ind_composites = v_ind_composites
    ),
    "Some variables in v_dep are not in data"
  )
  expect_error(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = c(v_ind_composites, "fred")
    ),
    "Some variables in v_ind_composites are not in data"
  )
})


test_that("v_ind, v_ind_composites, and v_dep are in colnames of R", {

  d1 <- d %>% mutate(fred = 1)


  expect_error(
    cond_maha(
      data = d1,
      R = R,
      v_dep = v_dep,
      v_ind = "fred",
      v_ind_composites = v_ind_composites
    ),
    "Some variables in v_ind are not in R: fred"
  )



  expect_error(
    cond_maha(
      data = d1,
      R = R,
      v_dep = c(v_dep, "fred"),
      v_ind = NULL,
      v_ind_composites = v_ind_composites
    ),
    "Some variables in v_dep are not in R: fred"
  )

  expect_error(
    cond_maha(
      data = d1,
      R = R,
      v_dep = v_dep,
      v_ind_composites = c(v_ind_composites, "fred")
    ),
    "Some variables in v_ind_composites are not in R: fred"
  )
})

test_that("v_ind and v_dep have no overlapping names", {
  expect_error(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind <- "X_1",
      v_ind_composites = v_ind_composites
    )
  )
  expect_error(cond_maha(
    data = d,
    R = R,
    v_dep = c(v_dep, "X"),
    v_ind_composites = v_ind_composites
  ))
  expect_error(cond_maha(
    data = d,
    R = R,
    v_dep = v_dep,
    v_ind_composites = c(v_ind_composites, "X_1")
  ))
})

test_that("Sample stats need at least 3 rows of data.", {
  expect_error(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites,
      use_sample_stats = TRUE
    )
  )
})

# Singularity ----

test_that("Measures not singular", {
  R_singular_dep <- R
  R_singular_dep["Y_1", "Y_2"] <- R_singular_dep["Y_2", "Y_1"] <- 1
  R_singular_ind <- R
  R_singular_ind["X_Composite", "Y_Composite"] <-
    R_singular_ind["Y_Composite", "X_Composite"] <- 1
  expect_error(
    cond_maha(
      data = d,
      R = R_singular_dep,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    )
  )
  expect_error(
    cond_maha(
      data = d,
      R = R_singular_ind,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    )
  )

  expect_true(is_singular(matrix(1,2,2)))
})

# Variable Metrics ----

test_that("Metric does not matter", {
  expect_equal(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    )$dCM,
    cond_maha(
      data = d %>% mutate_all(function(x) x * 15 + 100),
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites,
      mu = 100,
      sigma = 15
    )$dCM
  )
})

# Rounding ----
test_that("Propround", {
  expect_equal(
    proportion_round(
      c(-1,-.0001, 0,.000012,0.012,0.12,0.98,0.99,0.991,0.999,0.9991, 1, 1.001,2,NA)),
    c("-1.00","-0.00",".00",".000012",".01",".12",".98",".99",".991",".999",".9991","1.00","1.00","2.00", NA_character_))

  expect_equal(
    proportion_round(
      c(-1,-.00012, 0,.00001234,0.01234,0.1234,0.9111, 0.99,0.991111,0.999,0.9991111, 1, 1.001234,2,NA), digits = 3),
    c("-1.000","-0.000",".000",".000012",".012",".123",".911",".990",".991",".999",".9991","1.000","1.001","2.000", NA_character_))

  expect_equal(
    proportion2percentile(c(0.001,0.01,.5,.99,.992)),
    c(".1","1","50","99","99.2")
  )
})

test_that("Labeling",{
  expect_equal(
    p2label(
      pnorm(
        seq(60,140,10),
        mean = 100,
        sd = 15)),
    c("Extremely Low",
      "Very Low",
      "Low",
      "Low Average",
      "Average",
      "High Average",
      "High",
      "Very High",
      "Extremely High"))
})






data <- bind_rows(d,d %>% mutate(across(everything(), function(x) x - 1)))

cm <- cond_maha(
  data = data,
  R = R,
  v_dep = v_dep,
  # v_ind_composites = v_ind_composites
)
use_sample_stats <- FALSE
mu <- 0
sigma <- 1
library(ggplot2)
p_tail <- 0.40
plot(cm, p_tail = 0.05)

# m <- "
# Gc =~ 0.85 * Gc_1 + 0.68 * Gc_2 + 0.80 * Gc_3
# Gf =~ 0.80 * Gf_1 + 0.90 * Gf_2 + 0.80 * Gf_3
# Gs =~ 0.70 * Gs_1 + 0.80 * Gs_2 + 0.80 * Gs_3
# Read =~ 0.66 * Read_1 + 0.85 * Read_2 + 0.91 * Read_3
# Math =~ 0.40 * Math_1 + 0.90 * Math_2 + 0.70 * Math3
# Gc ~ 0.60 * Gf + 0.10 * Gs
# Gf ~ 0.50 * Gs
# Read ~ 0.40 * Gc + 0.10 * Gf
# Math ~ 0.20 * Gc + 0.30 * Gf + 0.10 * Gs
# "
# boot_dcm(model = m,
#          v_dep = c("Math", "Read"),
#          v_ind = c("Gc", "Gf", "Gs"),
#          sample_size = 100,
#          replications = 100)
# v_ind_composites <- NULL
