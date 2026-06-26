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

  out <- suppressWarnings(engine_brm(
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
  ))

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

test_that("brm measurement plan keeps plain leaf blocks as structural families", {
  plan <- make_brm_measurement_plan(
    measurement_blocks = list(alpha = c("a1", "a2"), beta = c("b1", "b2")),
    variables = c("a1", "a2", "b1", "b2")
  )

  expect_equal(names(plan$leaf_blocks), c("alpha", "beta"))
  expect_equal(names(plan$structural_families), c("alpha", "beta"))
  expect_equal(plan$structural_families$alpha$leaf_blocks, "alpha")
  expect_equal(plan$ignored_composites, character())
})

test_that("brm measurement plan preserves lower blocks and ignores covariate-only composites", {
  composite_list <- list(
    lower = list(
      ses = c("mom", "dad"),
      hardwork = c("h1", "h2"),
      ability = c("a1", "a2"),
      material = c("m1", "m2", "m3")
    ),
    higher = list(
      merit = c("hardwork", "ability")
    ),
    order = c("ses", "hardwork", "ability", "material", "merit")
  )

  plan <- make_brm_measurement_plan(
    measurement_blocks = composite_list,
    variables = c("h1", "h2", "a1", "a2", "m1", "m2", "m3")
  )

  expect_equal(names(plan$leaf_blocks), c("hardwork", "ability", "material"))
  expect_equal(names(plan$structural_families), c("merit", "material"))
  expect_equal(plan$structural_families$merit$leaf_blocks, c("hardwork", "ability"))
  expect_equal(plan$ignored_composites, "ses")
})

test_that("coherent gaussian epred sampler uses supplied residual donor rows", {
  epred_new <- rbind(c(4, 4, 4), c(8, 8, 8))
  epred_train <- matrix(4, nrow = 1, ncol = 3)
  y_train <- c(3, 4, 5)
  donor_rows <- c(3, 2, 1)
  epred_rows <- c(1, 2, 1)

  out <- with_seed(10, sample_coherent_gaussian_epred_column(
    epred_new = epred_new,
    epred_train = epred_train,
    y_train = y_train,
    template = y_train,
    donor_rows = donor_rows,
    epred_rows = epred_rows
  ))

  expect_equal(out, c(5, 5, 3))
})

test_that("measurement-aware brm engine uses block-score donor reconstruction", {
  skip_if_not_installed("brms")

  set.seed(42)
  latent <- rnorm(40)
  df <- data.frame(
    a1 = round(pmin(pmax(3 + latent + rnorm(40, sd = 0.2), 1), 5)),
    a2 = round(pmin(pmax(3 + latent + rnorm(40, sd = 0.2), 1), 5)),
    z = rnorm(40)
  )

  out <- engine_brm(
    data = df,
    variables = c("a1", "a2"),
    covariates = "z",
    measurement_blocks = list(alpha = c("a1", "a2")),
    between_block_shrink = 0,
    predictor_method = "full",
    predictor_k = 2,
    n_aug = 12,
    families = stats::setNames(rep(list(gaussian()), 2), c("a1", "a2")),
    iter = 20,
    warmup = 10,
    chains = 1,
    cores = 1,
    algorithm = "meanfield",
    seed = 2
  )

  expect_equal(nrow(out$data), 12)
  expect_equal(out$models$all$measurement_mode, "block_score_donor_brm")
  expect_equal(out$models$all$between_block_shrink, 0)
  expect_equal(names(out$models$all$leaf_blocks), "alpha")
  expect_equal(names(out$models$all$structural_families), "alpha")
  expect_equal(out$models$all$donor_settings$method, "nearest_leaf_score")
  expect_true(all(out$data$a1 %in% sort(unique(df$a1))))
  expect_true(all(out$data$a2 %in% sort(unique(df$a2))))
  expect_gt(stats::cor(out$data$a1, out$data$a2), 0.5)
})
