test_that("brm engine requires brms", {
  skip_if_not_installed("brms")
  expect_true(requireNamespace("brms", quietly = TRUE))
})

test_that("brm engine preserves continuous outcomes when using gaussian families", {
  skip_if_not_installed("brms")

  df <- data.frame(
    x = rnorm(30),
    y = rnorm(30),
    z = rnorm(30)
  )

  out <- engine_brm(
    data = df,
    variables = c("x", "y"),
    covariates = "z",
    n_aug = 10,
    iter = 20,
    warmup = 10,
    chains = 1,
    cores = 1,
    algorithm = "meanfield",
    seed = 1
  )

  expect_equal(nrow(out$data), 10)
  expect_type(out$data$x, "double")
  expect_type(out$data$y, "double")
})

test_that("measurement-aware brm role priors use expected scales", {
  skip_if_not_installed("brms")

  df <- data.frame(
    a2 = rnorm(20),
    b1 = rnorm(20),
    z = rnorm(20)
  )
  roles <- list(
    within_block_predictors = "a2",
    between_block_predictors = "b1",
    between_k = 1,
    target_block = "alpha"
  )

  out <- make_brm_role_priors(
    data = df,
    predictor_roles = roles,
    covariates = "z",
    base_sd = 10,
    shrink_sd = 2,
    between_block_shrink = 0.25
  )
  prior_text <- paste(capture.output(print(out$priors)), collapse = "\n")

  expect_equal(out$meta$prior_mode, "role_specific")
  expect_equal(out$meta$same_block_sd, 10)
  expect_equal(out$meta$between_block_sd, 0.5)
  expect_equal(out$meta$covariate_sd, 2)
  expect_match(prior_text, "normal\\(0, 10\\)")
  expect_match(prior_text, "normal\\(0, 0.5\\)")
  expect_match(prior_text, "normal\\(0, 2\\)")
})

test_that("measurement-aware brm prior builder falls back for unsafe terms", {
  skip_if_not_installed("brms")

  df <- data.frame(
    safe = rnorm(20),
    group = rep(c("a", "b"), 10)
  )
  roles <- list(
    within_block_predictors = "safe",
    between_block_predictors = character(),
    between_k = 0,
    target_block = "alpha"
  )

  out <- make_brm_role_priors(
    data = df,
    predictor_roles = roles,
    covariates = "group",
    base_sd = 10,
    shrink_sd = 2,
    between_block_shrink = 0.25
  )

  expect_equal(out$meta$prior_mode, "global_fallback")
  expect_equal(out$meta$fallback_terms, "group")
})

test_that("coherent gaussian epred sampler uses supplied residual donor rows", {
  epred_new <- matrix(4, nrow = 1, ncol = 3)
  epred_train <- matrix(4, nrow = 1, ncol = 3)
  y_train <- c(3, 4, 5)
  donor_rows <- c(3, 2, 1)

  out <- with_seed(10, sample_coherent_gaussian_epred_column(
    epred_new = epred_new,
    epred_train = epred_train,
    y_train = y_train,
    template = y_train,
    donor_rows = donor_rows
  ))

  expect_equal(out, c(5, 4, 3))
})

test_that("measurement-aware brm engine records predictor roles", {
  skip_if_not_installed("brms")

  set.seed(42)
  df <- data.frame(
    a1 = rnorm(24),
    a2 = rnorm(24),
    z = rnorm(24)
  )

  out <- engine_brm(
    data = df,
    variables = c("a1", "a2"),
    covariates = "z",
    measurement_blocks = list(alpha = c("a1", "a2")),
    between_block_shrink = 0,
    predictor_method = "full",
    predictor_k = 2,
    n_aug = 4,
    families = stats::setNames(rep(list(gaussian()), 2), c("a1", "a2")),
    iter = 20,
    warmup = 10,
    chains = 1,
    cores = 1,
    algorithm = "meanfield",
    seed = 2
  )

  expect_equal(nrow(out$data), 4)
  expect_equal(out$models$a1$measurement_mode, "block_predictor_brm")
  expect_equal(out$models$a1$between_block_shrink, 0)
  expect_equal(out$models$a1$predictor_roles$within_block_predictors, "a2")
  expect_equal(out$models$a1$predictor_roles$between_block_predictors, character())
  expect_equal(out$models$a1$predictor_roles$covariates, "z")
  expect_equal(out$models$a1$prior_roles$prior_mode, "role_specific")
})
