make_engine_fixture <- function() {
  data.frame(
    x1 = c(1, 2, 3, 4, 5, 2, 3, 4, 2, 5),
    x2 = c(5, 4, 3, 2, 1, 4, 3, 2, 4, 1),
    x3 = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
    x4 = c(1.2, 2.4, 3.1, 1.8, 2.7, 3.5, 2.2, 1.1, 2.0, 2.9),
    g1 = c("a", "a", "b", "b", "a", "a", "b", "b", "a", "b"),
    g2 = c(1, 1, 1, 1, 2, 2, 2, 2, 1, 2)
  )
}

test_that("sydni exposes engine_control in its public interface", {
  expect_true("engine_control" %in% names(formals(sydni)))
})

test_that("bart engine fails clearly when dbarts is unavailable", {
  if (requireNamespace("dbarts", quietly = TRUE)) {
    skip("dbarts is installed; unavailable-package test not applicable.")
  }

  dat <- make_engine_fixture()

  expect_error(
    engine_bart(
      data = dat,
      variables = c("x1", "x4"),
      covariates = c("g1", "g2"),
      n_aug = 4,
      seed = 1
    ),
    "dbarts"
  )
})

test_that("factor engine fails clearly when psych is unavailable", {
  if (requireNamespace("psych", quietly = TRUE)) {
    skip("psych is installed; unavailable-package test not applicable.")
  }

  dat <- make_engine_fixture()

  expect_error(
    engine_factor(
      data = dat,
      variables = c("x1", "x2", "x3"),
      covariates = c("g1", "g2"),
      n_aug = 4,
      seed = 1
    ),
    "psych"
  )
})

test_that("graphical engine fails clearly when huge is unavailable", {
  if (requireNamespace("huge", quietly = TRUE)) {
    skip("huge is installed; unavailable-package test not applicable.")
  }

  dat <- make_engine_fixture()

  expect_error(
    engine_graphical(
      data = dat,
      variables = c("x1", "x2", "x3", "x4"),
      covariates = c("g1", "g2"),
      n_aug = 4,
      seed = 1
    ),
    "huge"
  )
})

test_that("bart measurement plans preserve hierarchy and reject overlap", {
  composite_list <- list(
    lower = list(
      ses = c("mom", "dad"),
      material = c("m1", "m2"),
      luck = c("l1", "l2"),
      fate = c("f1", "f2")
    ),
    higher = list(
      chance = c("luck", "fate")
    ),
    order = c("ses", "material", "luck", "fate", "chance")
  )
  plan <- make_measurement_score_plan(
    measurement_blocks = composite_list,
    variables = c("m1", "m2", "l1", "l2", "f1", "f2")
  )

  expect_equal(names(plan$leaf_blocks), c("material", "luck", "fate"))
  expect_equal(plan$family_names, c("material", "chance"))
  expect_equal(plan$structural_families$chance$leaf_blocks, c("luck", "fate"))
  expect_equal(plan$ignored_composites, "ses")

  overlapping <- composite_list
  overlapping$higher$fortune <- c("luck", "material")
  expect_error(
    make_measurement_score_plan(
      measurement_blocks = overlapping,
      variables = c("m1", "m2", "l1", "l2", "f1", "f2")
    ),
    "only one final structural family"
  )
})

test_that("bart bundle calibration shrinks between families without changing reliability", {
  set.seed(110)
  n <- 300
  z1 <- rnorm(n)
  z2 <- 0.8 * z1 + rnorm(n, sd = 0.5)
  z3 <- -0.6 * z1 + rnorm(n, sd = 0.7)
  dat <- data.frame(
    a1 = z1 + rnorm(n, sd = 0.35),
    a2 = z1 + rnorm(n, sd = 0.35),
    b1 = z2 + rnorm(n, sd = 0.35),
    b2 = z2 + rnorm(n, sd = 0.35),
    b3 = z2 + rnorm(n, sd = 0.35),
    c1 = z3 + rnorm(n, sd = 0.35),
    c2 = z3 + rnorm(n, sd = 0.35)
  )
  plan <- make_measurement_score_plan(
    measurement_blocks = list(
      alpha = c("a1", "a2"),
      beta = c("b1", "b2", "b3"),
      gamma = c("c1", "c2")
    ),
    variables = names(dat)
  )
  target <- diag(3)
  dimnames(target) <- list(plan$family_names, plan$family_names)

  calibrated <- calibrate_measurement_structural_bundles(
    data = dat,
    measurement_plan = plan,
    target_cor = target,
    seed = 22,
    calibration_draws = 25
  )
  diagnostics <- attr(calibrated, "calibration")

  cronbach_alpha <- function(x) {
    x <- as.data.frame(x)
    p <- ncol(x)
    p / (p - 1) * (1 - sum(vapply(x, stats::var, numeric(1))) / stats::var(rowSums(x)))
  }
  row_keys <- function(x) {
    apply(x, 1, function(row) paste(format(row, digits = 17), collapse = "|"))
  }

  for (family_name in plan$family_names) {
    family_vars <- plan$structural_families[[family_name]]$variables
    expect_equal(
      sort(row_keys(dat[, family_vars, drop = FALSE])),
      sort(row_keys(calibrated[, family_vars, drop = FALSE]))
    )
    expect_equal(
      stats::cor(dat[, family_vars, drop = FALSE]),
      stats::cor(calibrated[, family_vars, drop = FALSE]),
      tolerance = 1e-12
    )
    expect_equal(
      cronbach_alpha(dat[, family_vars, drop = FALSE]),
      cronbach_alpha(calibrated[, family_vars, drop = FALSE]),
      tolerance = 1e-12
    )
  }

  expect_lte(diagnostics$error, diagnostics$error_raw)
  expect_lte(diagnostics$max_error, diagnostics$max_error_raw)
  expect_lt(diagnostics$max_error, 0.15)
  expect_lt(diagnostics$within_family_correlation_max_change, 1e-12)
  expect_equal(
    calibrated,
    calibrate_measurement_structural_bundles(
      data = dat,
      measurement_plan = plan,
      target_cor = target,
      seed = 22,
      calibration_draws = 25
    ),
    ignore_attr = TRUE
  )
})

test_that("bart bundle calibration approaches an observed correlation target", {
  set.seed(111)
  n <- 250
  z1 <- rnorm(n)
  z2 <- 0.7 * z1 + rnorm(n, sd = 0.7)
  observed <- data.frame(
    a1 = z1 + rnorm(n, sd = 0.3),
    a2 = z1 + rnorm(n, sd = 0.3),
    b1 = z2 + rnorm(n, sd = 0.3),
    b2 = z2 + rnorm(n, sd = 0.3)
  )
  generated <- observed
  generated[, c("b1", "b2")] <- generated[sample.int(n), c("b1", "b2")]
  plan <- make_measurement_score_plan(
    measurement_blocks = list(alpha = c("a1", "a2"), beta = c("b1", "b2")),
    variables = names(observed)
  )
  target <- estimate_brm_score_correlation(
    make_measurement_structural_score_frame(observed, plan)
  )
  calibrated <- calibrate_measurement_structural_bundles(
    data = generated,
    measurement_plan = plan,
    target_cor = target,
    seed = 23,
    calibration_draws = 25
  )
  diagnostics <- attr(calibrated, "calibration")

  expect_lt(diagnostics$error, diagnostics$error_raw)
  expect_lt(diagnostics$max_error, 0.15)
})

test_that("bart engine returns models and data when dbarts is installed", {
  skip_if_not_installed("dbarts")

  dat <- make_engine_fixture()

  seen_warning <- FALSE
  res <- withCallingHandlers(
    engine_bart(
      data = dat,
      variables = c("x1", "x2", "x3", "x4"),
      covariates = c("g1", "g2"),
      n_aug = 6,
      predictor_k = 2,
      engine_control = list(num_trees = 10, num_burn = 20, num_draws = 20),
      seed = 1
    ),
    warning = function(w) {
      if (grepl("ordinal-like", conditionMessage(w), fixed = TRUE)) {
        seen_warning <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )

  expect_true(seen_warning)
  expect_named(res, c("models", "data"))
  expect_equal(nrow(res$data), 6)
  expect_equal(names(res$data), c("x1", "x2", "x3", "x4", "g1", "g2"))
})

test_that("bart ordinal posterior predictive draws do not collapse to one category", {
  skip_if_not_installed("dbarts")

  set.seed(31)
  n <- 120
  latent <- rnorm(n)
  dat <- data.frame(
    x1 = as.integer(cut(latent + rnorm(n, sd = 0.8), 5, labels = FALSE)),
    x2 = as.integer(cut(latent + rnorm(n, sd = 0.8), 5, labels = FALSE)),
    x3 = as.integer(cut(rnorm(n), 5, labels = FALSE))
  )

  res <- suppressWarnings(engine_bart(
    data = dat,
    variables = names(dat),
    predictor_k = 2,
    n_aug = 80,
    engine_control = list(num_trees = 5, num_burn = 5, num_draws = 5),
    seed = 31
  ))

  expect_gt(length(unique(res$data$x1)), 1)
  expect_true(all(res$data$x1 %in% sort(unique(dat$x1))))
})

test_that("measurement-aware bart engine records predictor roles", {
  skip_if_not_installed("dbarts")

  set.seed(41)
  n <- 90
  latent_a <- rnorm(n)
  latent_b <- 0.5 * latent_a + rnorm(n, sd = 0.9)
  dat <- data.frame(
    a1 = as.integer(cut(latent_a + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    a2 = as.integer(cut(latent_a + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    a3 = as.integer(cut(latent_a + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    b1 = as.integer(cut(latent_b + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    b2 = as.integer(cut(latent_b + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    b3 = as.integer(cut(latent_b + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    g = rep(c("x", "y"), length.out = n)
  )

  res <- suppressWarnings(engine_bart(
    data = dat,
    variables = c("a1", "a2", "a3", "b1", "b2", "b3"),
    covariates = "g",
    measurement_blocks = list(alpha = c("a1", "a2", "a3"), beta = c("b1", "b2", "b3")),
    between_block_shrink = 0.2,
    predictor_method = "full",
    predictor_k = 3,
    n_aug = 20,
    engine_control = list(num_trees = 5, num_burn = 5, num_draws = 5),
    seed = 6
  ))

  expect_equal(nrow(res$data), 20)
  expect_true(all(c("a1", "a2", "a3", "b1", "b2", "b3", "g") %in% names(res$data)))
  expect_true(all(res$data$a1 %in% sort(unique(dat$a1))))
  expect_equal(res$models$a1$measurement_mode, "block_bundle_calibrated_bart")
  expect_equal(res$models$a1$between_block_shrink, 0.2)
  expect_equal(names(res$models$a1$measurement_blocks), c("alpha", "beta"))
  expect_equal(sort(res$models$a1$predictor_roles$within_block_predictors), c("a2", "a3"))
  expect_equal(length(res$models$a1$predictor_roles$between_block_predictors), 0)
  expect_equal(res$models$a1$predictor_roles$covariates, "g")
  expect_equal(res$models$all$measurement_mode, "block_bundle_calibrated_bart")
  expect_equal(res$models$all$score_calibration_mode, "iman_conover_family_bundle")
  expect_equal(res$models$all$score_calibration_draws, 25)
  expect_lte(
    res$models$all$structural_correlation_error_calibrated,
    res$models$all$structural_correlation_error_raw
  )
  expect_lte(
    res$models$all$structural_correlation_max_error_calibrated,
    res$models$all$structural_correlation_max_error_raw
  )
  expect_lt(res$models$all$within_family_correlation_max_change, 1e-12)
  expect_equal(dim(res$models$all$structural_covariate_correlation_generated), c(2, 2))
})

test_that("measurement-aware bart between-block predictors follow shrinkage", {
  skip_if_not_installed("dbarts")

  dat <- data.frame(
    a1 = rnorm(50),
    a2 = rnorm(50),
    a3 = rnorm(50),
    b1 = rnorm(50),
    b2 = rnorm(50),
    b3 = rnorm(50)
  )

  low <- engine_bart(
    data = dat,
    variables = names(dat),
    measurement_blocks = list(alpha = c("a1", "a2", "a3"), beta = c("b1", "b2", "b3")),
    between_block_shrink = 0,
    predictor_method = "full",
    predictor_k = 3,
    n_aug = 10,
    engine_control = list(num_trees = 5, num_burn = 5, num_draws = 5),
    seed = 7
  )
  high <- engine_bart(
    data = dat,
    variables = names(dat),
    measurement_blocks = list(alpha = c("a1", "a2", "a3"), beta = c("b1", "b2", "b3")),
    between_block_shrink = 1,
    predictor_method = "full",
    predictor_k = 3,
    n_aug = 10,
    engine_control = list(num_trees = 5, num_burn = 5, num_draws = 5),
    seed = 7
  )

  expect_equal(length(low$models$a1$predictor_roles$between_block_predictors), 0)
  expect_equal(length(high$models$a1$predictor_roles$between_block_predictors), 3)
})

test_that("measurement-aware bart accepts cscore-style higher-order blocks", {
  skip_if_not_installed("dbarts")

  dat <- data.frame(
    gs1 = rnorm(50),
    gs2 = rnorm(50),
    gs3 = rnorm(50),
    gs4 = rnorm(50),
    e1 = rnorm(50),
    e2 = rnorm(50),
    mom = rnorm(50),
    dad = rnorm(50)
  )
  composite_list <- list(
    lower = list(
      ses = c("mom", "dad"),
      consistency_interest = c("gs2", "gs3"),
      perseverance_effort = c("gs1", "gs4"),
      extraversion = c("e1", "e2")
    ),
    higher = list(
      grit = c("consistency_interest", "perseverance_effort")
    ),
    order = c("ses", "consistency_interest", "perseverance_effort", "extraversion", "grit")
  )

  res <- engine_bart(
    data = dat,
    variables = setdiff(names(dat), c("mom", "dad")),
    covariates = c("mom", "dad"),
    measurement_blocks = composite_list,
    n_aug = 10,
    predictor_k = 2,
    engine_control = list(num_trees = 5, num_burn = 5, num_draws = 5),
    seed = 8
  )

  expect_equal(names(res$models$gs1$measurement_blocks), c("extraversion", "grit"))
  expect_equal(sort(res$models$gs1$measurement_blocks$grit), sort(c("gs1", "gs2", "gs3", "gs4")))
  expect_equal(res$models$all$ignored_composites, "ses")
  expect_true(all(res$data$mom %in% dat$mom))
  expect_true(all(res$data$dad %in% dat$dad))
})

test_that("measurement-aware bart keeps same-block predictors when predictor_k is zero", {
  skip_if_not_installed("dbarts")

  dat <- data.frame(
    a1 = rnorm(40),
    a2 = rnorm(40),
    b1 = rnorm(40),
    b2 = rnorm(40)
  )

  res <- engine_bart(
    data = dat,
    variables = names(dat),
    measurement_blocks = list(alpha = c("a1", "a2"), beta = c("b1", "b2")),
    predictor_k = 0,
    n_aug = 8,
    engine_control = list(num_trees = 5, num_burn = 5, num_draws = 5),
    seed = 9
  )

  expect_equal(res$models$a1$predictor_roles$within_block_predictors, "a2")
  expect_equal(res$models$a1$predictor_roles$between_block_predictors, character())
})

test_that("factor engine subgroups by covariate combinations when psych is installed", {
  skip_if_not_installed("psych")

  dat <- make_engine_fixture()

  res <- engine_factor(
    data = dat,
    variables = c("x1", "x2", "x3"),
    covariates = c("g1", "g2"),
    n_aug = 8,
    engine_control = list(n_factors = 1),
    seed = 1
  )

  expect_named(res, c("models", "data"))
  expect_equal(sort(unique(paste(res$data$g1, res$data$g2))), sort(unique(paste(dat$g1, dat$g2))))
})

test_that("measurement-aware factor engine preserves block metadata and item support", {
  skip_if_not_installed("psych")

  set.seed(12)
  n <- 120
  latent_a <- rnorm(n)
  latent_b <- 0.6 * latent_a + rnorm(n, sd = 0.8)
  dat <- data.frame(
    a1 = as.integer(cut(latent_a + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    a2 = as.integer(cut(latent_a + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    a3 = as.integer(cut(latent_a + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    b1 = as.integer(cut(latent_b + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    b2 = as.integer(cut(latent_b + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    b3 = as.integer(cut(latent_b + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    g = rep(c("x", "y"), each = n / 2)
  )

  res <- engine_factor(
    data = dat,
    variables = c("a1", "a2", "a3", "b1", "b2", "b3"),
    covariates = "g",
    measurement_blocks = list(alpha = c("a1", "a2", "a3"), beta = c("b1", "b2", "b3")),
    between_block_shrink = 0.2,
    n_aug = 40,
    seed = 4
  )

  expect_equal(nrow(res$data), 40)
  expect_true(all(c("a1", "a2", "a3", "b1", "b2", "b3", "g") %in% names(res$data)))
  expect_true(all(res$data$a1 %in% sort(unique(dat$a1))))
  expect_equal(res$models$x$measurement_mode, "block_one_factor_calibrated")
  expect_equal(res$models$x$between_block_shrink, 0.2)
  expect_equal(names(res$models$x$measurement_blocks), c("alpha", "beta"))
  expect_lt(
    abs(res$models$x$block_correlation_global["alpha", "beta"]),
    abs(res$models$x$block_correlation_global_raw["alpha", "beta"])
  )
  expect_equal(res$models$all$measurement_mode, "block_one_factor_calibrated")
  expect_equal(res$models$all$score_calibration_mode, "iman_conover_family_bundle")
  expect_equal(res$models$all$score_calibration_draws, 25)
  expect_lte(
    res$models$all$structural_correlation_error_calibrated,
    res$models$all$structural_correlation_error_raw
  )
  expect_lte(
    res$models$all$structural_correlation_max_error_calibrated,
    res$models$all$structural_correlation_max_error_raw
  )
  expect_lt(res$models$all$within_family_correlation_max_change, 1e-12)
  expect_equal(dim(res$models$all$structural_covariate_correlation_generated), c(2, 2))
})

test_that("measurement-aware factor engine accepts cscore-style higher-order blocks", {
  skip_if_not_installed("psych")

  set.seed(22)
  n <- 100
  latent_grit <- rnorm(n)
  latent_e <- 0.3 * latent_grit + rnorm(n, sd = 0.9)
  dat <- data.frame(
    gs1 = as.integer(cut(latent_grit + rnorm(n, sd = 0.4), 5, labels = FALSE)),
    gs2 = as.integer(cut(latent_grit + rnorm(n, sd = 0.4), 5, labels = FALSE)),
    gs3 = as.integer(cut(latent_grit + rnorm(n, sd = 0.4), 5, labels = FALSE)),
    gs4 = as.integer(cut(latent_grit + rnorm(n, sd = 0.4), 5, labels = FALSE)),
    e1 = as.integer(cut(latent_e + rnorm(n, sd = 0.4), 5, labels = FALSE)),
    e2 = as.integer(cut(latent_e + rnorm(n, sd = 0.4), 5, labels = FALSE)),
    mom = rep(c(1, 2), each = n / 2),
    dad = rep(rep(c(1, 2), each = n / 4), 2)
  )
  composite_list <- list(
    lower = list(
      ses = c("mom", "dad"),
      consistency_interest = c("gs2", "gs3"),
      perseverance_effort = c("gs1", "gs4"),
      extraversion = c("e1", "e2")
    ),
    higher = list(
      grit = c("consistency_interest", "perseverance_effort")
    ),
    order = c("ses", "consistency_interest", "perseverance_effort", "extraversion", "grit")
  )

  res <- engine_factor(
    data = dat,
    variables = setdiff(names(dat), c("mom", "dad")),
    covariates = c("mom", "dad"),
    measurement_blocks = composite_list,
    n_aug = 20,
    seed = 8
  )

  expect_equal(names(res$models$all$measurement_blocks), c("extraversion", "grit"))
  expect_equal(sort(res$models$all$measurement_blocks$grit), sort(c("gs1", "gs2", "gs3", "gs4")))
  expect_equal(res$models$all$ignored_composites, "ses")
  expect_equal(res$models$all$measurement_mode, "block_one_factor_calibrated")
  expect_true(all(res$data$mom %in% dat$mom))
  expect_true(all(res$data$dad %in% dat$dad))
})

test_that("graphical engine subgroups by covariate combinations when huge is installed", {
  skip_if_not_installed("huge")

  dat <- make_engine_fixture()

  res <- engine_graphical(
    data = dat,
    variables = c("x1", "x2", "x3", "x4"),
    covariates = c("g1", "g2"),
    n_aug = 8,
    seed = 1
  )

  expect_named(res, c("models", "data"))
  expect_equal(sort(unique(paste(res$data$g1, res$data$g2))), sort(unique(paste(dat$g1, dat$g2))))
})

test_that("measurement-aware graphical engine preserves block metadata and item support", {
  skip_if_not_installed("huge")

  set.seed(31)
  n <- 120
  latent_a <- rnorm(n)
  latent_b <- 0.6 * latent_a + rnorm(n, sd = 0.8)
  dat <- data.frame(
    a1 = as.integer(cut(latent_a + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    a2 = as.integer(cut(latent_a + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    a3 = as.integer(cut(latent_a + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    b1 = as.integer(cut(latent_b + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    b2 = as.integer(cut(latent_b + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    b3 = as.integer(cut(latent_b + rnorm(n, sd = 0.3), 5, labels = FALSE)),
    g = rep(c("x", "y"), each = n / 2)
  )

  res <- engine_graphical(
    data = dat,
    variables = c("a1", "a2", "a3", "b1", "b2", "b3"),
    covariates = "g",
    measurement_blocks = list(alpha = c("a1", "a2", "a3"), beta = c("b1", "b2", "b3")),
    between_block_shrink = 0.2,
    n_aug = 40,
    engine_control = list(transform = "scale", selection = "fixed_lambda", lambda_fixed = 0.01),
    seed = 5
  )

  expect_equal(nrow(res$data), 40)
  expect_true(all(c("a1", "a2", "a3", "b1", "b2", "b3", "g") %in% names(res$data)))
  expect_true(all(res$data$a1 %in% sort(unique(dat$a1))))
  expect_equal(res$models$x$measurement_mode, "block_graphical_calibrated")
  expect_equal(res$models$x$between_block_shrink, 0.2)
  expect_equal(names(res$models$x$measurement_blocks), c("alpha", "beta"))
  expect_gt(
    abs(res$models$x$correlation["a1", "a2"]),
    abs(res$models$x$correlation["a1", "b1"])
  )
  expect_lt(
    abs(res$models$x$correlation_measurement_adjusted["a1", "b1"]),
    abs(res$models$x$correlation_raw["a1", "b1"])
  )
  expect_equal(res$models$all$measurement_mode, "block_graphical_calibrated")
  expect_equal(res$models$all$score_calibration_mode, "iman_conover_family_bundle")
  expect_equal(res$models$all$score_calibration_draws, 25)
  expect_lte(
    res$models$all$structural_correlation_error_calibrated,
    res$models$all$structural_correlation_error_raw
  )
  expect_lte(
    res$models$all$structural_correlation_max_error_calibrated,
    res$models$all$structural_correlation_max_error_raw
  )
  expect_lt(res$models$all$within_family_correlation_max_change, 1e-12)
  expect_equal(dim(res$models$all$structural_covariate_correlation_generated), c(2, 2))
})

test_that("measurement-aware graphical engine accepts cscore-style higher-order blocks", {
  skip_if_not_installed("huge")

  set.seed(32)
  n <- 100
  latent_grit <- rnorm(n)
  latent_e <- 0.3 * latent_grit + rnorm(n, sd = 0.9)
  dat <- data.frame(
    gs1 = as.integer(cut(latent_grit + rnorm(n, sd = 0.4), 5, labels = FALSE)),
    gs2 = as.integer(cut(latent_grit + rnorm(n, sd = 0.4), 5, labels = FALSE)),
    gs3 = as.integer(cut(latent_grit + rnorm(n, sd = 0.4), 5, labels = FALSE)),
    gs4 = as.integer(cut(latent_grit + rnorm(n, sd = 0.4), 5, labels = FALSE)),
    e1 = as.integer(cut(latent_e + rnorm(n, sd = 0.4), 5, labels = FALSE)),
    e2 = as.integer(cut(latent_e + rnorm(n, sd = 0.4), 5, labels = FALSE)),
    mom = rep(c(1, 2), each = n / 2),
    dad = rep(rep(c(1, 2), each = n / 4), 2)
  )
  composite_list <- list(
    lower = list(
      ses = c("mom", "dad"),
      consistency_interest = c("gs2", "gs3"),
      perseverance_effort = c("gs1", "gs4"),
      extraversion = c("e1", "e2")
    ),
    higher = list(
      grit = c("consistency_interest", "perseverance_effort")
    ),
    order = c("ses", "consistency_interest", "perseverance_effort", "extraversion", "grit")
  )

  res <- engine_graphical(
    data = dat,
    variables = setdiff(names(dat), c("mom", "dad")),
    covariates = c("mom", "dad"),
    measurement_blocks = composite_list,
    n_aug = 20,
    engine_control = list(transform = "scale", selection = "fixed_lambda", lambda_fixed = 0.01),
    seed = 9
  )

  expect_equal(names(res$models$all$measurement_blocks), c("extraversion", "grit"))
  expect_equal(sort(res$models$all$measurement_blocks$grit), sort(c("gs1", "gs2", "gs3", "gs4")))
  expect_equal(res$models$all$measurement_mode, "block_graphical_calibrated")
  expect_equal(res$models$all$ignored_composites, "ses")
  expect_true(all(res$data$mom %in% dat$mom))
  expect_true(all(res$data$dad %in% dat$dad))
})

test_that("measurement-aware graphical fallback still records measurement metadata", {
  skip_if_not_installed("huge")

  dat <- data.frame(
    a1 = rep(1:5, 10),
    a2 = rep(3, 50)
  )

  res <- engine_graphical(
    data = dat,
    variables = c("a1", "a2"),
    measurement_blocks = list(alpha = c("a1", "a2")),
    n_aug = 10,
    engine_control = list(transform = "scale"),
    seed = 10
  )

  expect_true(res$models$all$fallback)
  expect_equal(res$models$all$measurement_mode, "block_graphical_calibrated")
  expect_equal(names(res$models$all$measurement_blocks), "alpha")
})
