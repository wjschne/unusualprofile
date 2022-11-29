# code to prepare `DATASET` dataset goes here

library(simstandard)
# lavaan model with three indicators of a latent variable
model <- "
X =~ 0.7 * X_1 + 0.5 * X_2 + 0.8 * X_3
Y =~ 0.8 * Y_1 + 0.7 * Y_2 + 0.9 * Y_3
Y ~ 0.6 * X
"

# Ensure that random data will be the same as in the example.
set.seed(281)

vnames <- c(
  "X_1", "X_2", "X_3",
  "Y_1", "Y_2", "Y_3",
  "X", "Y")

# Randomly generated case
d_example <- sim_standardized(
  model,
  n = 1,
  observed = TRUE,
  latent = TRUE,
  errors = FALSE,
  composites = FALSE
)

# Model-implied correlation matrix
R_example <- get_model_implied_correlations(model, latent = TRUE)

usethis::use_data(R_example, d_example, overwrite = TRUE)
