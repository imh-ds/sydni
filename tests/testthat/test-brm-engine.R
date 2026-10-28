test_that("brm engine fails fast when brms is unavailable", {
  isolated_engine <- engine_brm
  isolated_env <- new.env(parent = environment(engine_brm))
  isolated_env$requireNamespace <- function(package, quietly = FALSE) {
    if (identical(package, "brms")) return(FALSE)
    base::requireNamespace(package, quietly = quietly)
  }
  environment(isolated_engine) <- isolated_env
  testthat::local_mocked_bindings(engine_brm = isolated_engine, .package = "sydni")

  df <- data.frame(x = 1:5, y = 6:10)

  expect_error(
    isolated_engine(
      data = df,
      variables = c("x", "y"),
      n_aug = 3
    ),
    "The 'brms' package must be installed"
  )

  expect_error(
    sydni(
      data = df,
      variables = c("x", "y"),
      engine = "brm",
      n_aug = 3
    ),
    "The 'brms' package must be installed"
  )
})

test_that("empirical BRM intercept priors skip factor outcomes", {
  expect_null(estimate_intercept_prior(factor(c("no", "yes", "no"))))
  expect_null(estimate_intercept_prior(ordered(c("low", "high", "low"))))
  expect_null(estimate_intercept_prior(c(0L, 1L, 1L), family = stats::poisson()))
})

test_that("package metadata keeps brms optional", {
  desc <- readLines(testthat::test_path("..", "..", "DESCRIPTION"))
  namespace <- readLines(testthat::test_path("..", "..", "NAMESPACE"))

  expect_false(any(grepl("^\\s*brms\\s*,?$", desc[seq(which(grepl("^Imports:$", desc)) + 1, which(grepl("^Suggests:$", desc)) - 1)])))
  expect_true(any(grepl("^\\s*brms\\s*,?$", desc[seq(which(grepl("^Suggests:$", desc)) + 1, which(grepl("^Remotes:$", desc)) - 1)])))
  expect_false(any(grepl("^importFrom\\(brms,", namespace)))
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

test_that("brm structural residual donors are independent by family", {
  donors <- make_brm_structural_residual_donors(
    n_obs = 25,
    n_aug = 12,
    family_names = c("alpha", "beta", "gamma"),
    seed = 10
  )

  expect_equal(names(donors), c("alpha", "beta", "gamma"))
  expect_length(donors$alpha, 12)
  expect_true(all(donors$alpha >= 1 & donors$alpha <= 25))
  expect_false(identical(donors$alpha, donors$beta))
  expect_false(identical(donors$beta, donors$gamma))
})

test_that("brm structural score calibration preserves marginals and shrinks correlations", {
  set.seed(11)
  n <- 800
  x <- rnorm(n)
  raw_scores <- data.frame(
    alpha = x,
    beta = x + rnorm(n, sd = 0.1),
    gamma = rnorm(n)
  )
  target_cor <- diag(3)
  rownames(target_cor) <- colnames(target_cor) <- names(raw_scores)

  calibrated <- calibrate_brm_structural_score_correlations(
    score_data = raw_scores,
    target_cor = target_cor,
    seed = 12,
    calibration_draws = 25L
  )
  calibration <- attr(calibrated, "calibration")

  expect_equal(sort(calibrated$alpha), sort(raw_scores$alpha))
  expect_equal(sort(calibrated$beta), sort(raw_scores$beta))
  expect_equal(sort(calibrated$gamma), sort(raw_scores$gamma))
  expect_lt(calibration$error, calibration$raw_error)
  expect_lt(calibration$max_error, calibration$raw_max_error)
  expect_lt(
    abs(stats::cor(calibrated$alpha, calibrated$beta)),
    abs(stats::cor(raw_scores$alpha, raw_scores$beta))
  )
  expect_lt(abs(stats::cor(calibrated$alpha, calibrated$beta)), 0.15)
})

test_that("brm structural score calibration can target observed score correlations", {
  set.seed(13)
  n <- 1000
  target_cor <- matrix(
    c(
      1.00, 0.55, -0.25,
      0.55, 1.00, 0.20,
      -0.25, 0.20, 1.00
    ),
    nrow = 3,
    byrow = TRUE
  )
  colnames(target_cor) <- rownames(target_cor) <- c("alpha", "beta", "gamma")
  raw_scores <- data.frame(
    alpha = rt(n, df = 6),
    beta = rnorm(n),
    gamma = runif(n, min = -2, max = 2)
  )

  calibrated <- calibrate_brm_structural_score_correlations(
    score_data = raw_scores,
    target_cor = target_cor,
    seed = 14,
    calibration_draws = 25L
  )
  calibrated_cor <- stats::cor(calibrated)
  calibration <- attr(calibrated, "calibration")

  expect_equal(sort(calibrated$gamma), sort(raw_scores$gamma))
  expect_equal(calibration$draws, 25L)
  expect_lt(abs(calibrated_cor["alpha", "beta"] - target_cor["alpha", "beta"]), 0.15)
  expect_lt(abs(calibrated_cor["alpha", "gamma"] - target_cor["alpha", "gamma"]), 0.15)
  expect_lt(abs(calibrated_cor["beta", "gamma"] - target_cor["beta", "gamma"]), 0.15)
})

test_that("brm structural score calibration selects the best candidate draw", {
  set.seed(15)
  n <- 600
  raw_scores <- data.frame(
    alpha = rnorm(n),
    beta = rexp(n),
    gamma = stats::rt(n, df = 5)
  )
  target_cor <- matrix(
    c(
      1.0, 0.4, 0.2,
      0.4, 1.0, -0.3,
      0.2, -0.3, 1.0
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(target_cor) <- colnames(target_cor) <- names(raw_scores)

  one_draw <- calibrate_brm_structural_score_correlations(
    score_data = raw_scores,
    target_cor = target_cor,
    seed = 16,
    calibration_draws = 1L
  )
  many_draws <- calibrate_brm_structural_score_correlations(
    score_data = raw_scores,
    target_cor = target_cor,
    seed = 16,
    calibration_draws = 25L
  )
  one_meta <- attr(one_draw, "calibration")
  many_meta <- attr(many_draws, "calibration")

  expect_true(many_meta$max_error <= one_meta$max_error + 1e-12)
  expect_equal(many_meta$draws, 25L)
  expect_true(many_meta$selected_draw %in% seq_len(25L))
})

test_that("two-item brm donor reconstruction uses rank-matched full item pairs", {
  data <- data.frame(
    a1 = c(1, 2, 4, 5),
    a2 = c(1, 2, 4, 5),
    b1 = c(1, 3, 5, 7),
    b2 = c(2, 4, 6, 8),
    b3 = c(3, 5, 7, 9)
  )
  observed_leaf_scores <- data.frame(
    alpha = c(-2, -1, 1, 2),
    beta = c(-2, -1, 1, 2)
  )
  synthetic_leaf_scores <- data.frame(
    alpha = c(2, -2, 1, -1),
    beta = c(2, -2, 1, -1)
  )

  out <- reconstruct_brm_items_from_leaf_scores(
    synthetic_leaf_scores = synthetic_leaf_scores,
    observed_leaf_scores = observed_leaf_scores,
    data = data,
    leaf_blocks = list(alpha = c("a1", "a2"), beta = c("b1", "b2", "b3")),
    variables = names(data),
    seed = 1,
    donor_k = 1
  )
  donors <- attr(out, "donor_rows")

  expect_equal(donors$alpha, c(4L, 1L, 3L, 2L))
  expect_equal(out$a1, data$a1[donors$alpha])
  expect_equal(out$a2, data$a2[donors$alpha])
  expect_equal(out$b1, data$b1[donors$beta])
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
  expect_equal(out$models$all$donor_settings$two_item_method, "rank_quantile_leaf_score")
  expect_equal(out$models$all$residual_donor_mode, "independent_by_structural_family")
  expect_equal(out$models$all$score_calibration_mode, "iman_conover_rank_target")
  expect_equal(out$models$all$score_calibration_draws, 25L)
  expect_true(is.matrix(out$models$all$structural_correlation_target))
  expect_equal(dim(out$models$all$structural_correlation_target), c(1L, 1L))
  expect_equal(out$models$all$structural_correlation_error_raw, 0)
  expect_equal(out$models$all$structural_correlation_error_calibrated, 0)
  expect_equal(out$models$all$structural_correlation_max_error_raw, 0)
  expect_equal(out$models$all$structural_correlation_max_error_calibrated, 0)
  expect_true(all(out$data$a1 %in% sort(unique(df$a1))))
  expect_true(all(out$data$a2 %in% sort(unique(df$a2))))
  expect_gt(stats::cor(out$data$a1, out$data$a2), 0.5)
})
