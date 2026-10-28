test_that("copula engine preserves requested row count and covariates", {
  df <- data.frame(
    x = rep(1:5, 20),
    y = rep(5:1, 20),
    g = rep(c("A", "B"), each = 50)
  )

  out <- engine_copula(
    data = df,
    variables = c("x", "y"),
    covariates = "g",
    n_aug = 40,
    seed = 42
  )

  expect_equal(nrow(out$data), 40)
  expect_true(all(c("x", "y", "g") %in% names(out$data)))
  expect_true(all(out$data$x >= min(df$x) & out$data$x <= max(df$x)))
})

test_that("copula engine keeps integer-like synthesized variables integer-like", {
  df <- data.frame(
    x = rep(1:5, 20),
    y = rep(2:6, 20)
  )

  out <- engine_copula(
    data = df,
    variables = c("x", "y"),
    n_aug = 30,
    seed = 7
  )

  expect_true(is.integer(out$data$x) || all(out$data$x == round(out$data$x)))
  expect_true(is.integer(out$data$y) || all(out$data$y == round(out$data$y)))
  expect_null(out$models$all$measurement_mode)
  expect_null(out$models$all$score_calibration_mode)
})

test_that("explicit shrinkage without blocks creates singleton families", {
  set.seed(40)
  n <- 180
  x <- rnorm(n)
  y <- 0.7 * x + rnorm(n, sd = 0.7)
  z <- -0.4 * x + rnorm(n, sd = 0.9)
  df <- data.frame(x = x, y = y, z = z)

  out <- engine_copula(
    data = df,
    variables = names(df),
    between_block_shrink = 0.35,
    n_aug = 120,
    seed = 41
  )
  meta <- out$models$all

  expect_equal(meta$measurement_mode, "singleton_copula_calibrated")
  expect_equal(meta$measurement_blocks_source, "implicit_singletons")
  expect_equal(meta$shrinkage_scope, "all_synthesized_variables")
  expect_equal(names(meta$measurement_blocks), names(df))
  expect_true(all(vapply(meta$measurement_blocks, length, integer(1)) == 1L))
  expect_equal(
    meta$structural_correlation_target[row(meta$structural_correlation_target) != col(meta$structural_correlation_target)],
    0.35 * meta$structural_correlation_observed[row(meta$structural_correlation_observed) != col(meta$structural_correlation_observed)],
    tolerance = 1e-8
  )
  expect_lte(
    meta$structural_correlation_error_calibrated,
    meta$structural_correlation_error_raw
  )
  expect_lte(
    meta$structural_correlation_max_error_calibrated,
    meta$structural_correlation_max_error_raw
  )
  expect_equal(nrow(out$data), 120)
  expect_true(all(out$data$x >= min(df$x) & out$data$x <= max(df$x)))
})

test_that("singleton copula supports zero and observed-correlation targets", {
  set.seed(42)
  n <- 140
  x <- rnorm(n)
  df <- data.frame(
    x = x,
    y = 0.6 * x + rnorm(n, sd = 0.8),
    z = -0.5 * x + rnorm(n, sd = 0.8)
  )

  zero <- engine_copula(
    data = df,
    variables = names(df),
    between_block_shrink = 0,
    n_aug = 80,
    seed = 43
  )
  one <- engine_copula(
    data = df,
    variables = names(df),
    between_block_shrink = 1,
    n_aug = 80,
    seed = 43
  )

  zero_target <- zero$models$all$structural_correlation_target
  expect_lt(max(abs(zero_target[row(zero_target) != col(zero_target)])), 1e-8)
  expect_equal(
    one$models$all$structural_correlation_target,
    one$models$all$structural_correlation_observed,
    tolerance = 1e-8
  )
})

test_that("partial copula blocks add standalone variables as singleton families", {
  set.seed(44)
  n <- 160
  latent <- rnorm(n)
  df <- data.frame(
    a1 = latent + rnorm(n, sd = 0.3),
    a2 = latent + rnorm(n, sd = 0.3),
    x = 0.5 * latent + rnorm(n, sd = 0.8),
    y = rnorm(n)
  )

  out <- engine_copula(
    data = df,
    variables = names(df),
    measurement_blocks = list(alpha = c("a1", "a2")),
    between_block_shrink = 0.25,
    n_aug = 100,
    seed = 45
  )
  meta <- out$models$all

  expect_equal(names(meta$measurement_blocks), c("alpha", "x", "y"))
  expect_equal(meta$measurement_blocks$x, "x")
  expect_equal(meta$measurement_blocks$y, "y")
  expect_equal(meta$measurement_mode, "block_copula_calibrated")
  expect_equal(meta$measurement_blocks_source, "user_supplied")
  expect_equal(meta$shrinkage_scope, "between_structural_families")
  expect_lt(meta$within_family_correlation_max_change, 1e-12)
})

test_that("implicit singleton shrinkage validates its explicit multiplier", {
  df <- data.frame(x = rnorm(30), y = rnorm(30))

  expect_error(
    engine_copula(
      data = df,
      variables = names(df),
      between_block_shrink = 1.2,
      n_aug = 10,
      seed = 46
    ),
    "must be a single numeric value in \\[0, 1\\]"
  )
})

test_that("stabilized copula simulation handles nearly singular correlations", {
  sigma <- matrix(c(
    1, 0.9999, 0.9998,
    0.9999, 1, 0.9997,
    0.9998, 0.9997, 1
  ), nrow = 3, byrow = TRUE)

  draws <- simulate_gaussian_copula(10, sigma)

  expect_equal(dim(draws), c(10, 3))
  expect_true(all(draws >= 0 & draws <= 1))
})

test_that("copula engine applies measurement-aware block shrinkage", {
  set.seed(1)
  latent_a <- rnorm(120)
  latent_b <- 0.7 * latent_a + rnorm(120, sd = 0.7)
  df <- data.frame(
    a1 = as.integer(cut(latent_a + rnorm(120, sd = 0.3), 5, labels = FALSE)),
    a2 = as.integer(cut(latent_a + rnorm(120, sd = 0.3), 5, labels = FALSE)),
    b1 = as.integer(cut(latent_b + rnorm(120, sd = 0.3), 5, labels = FALSE)),
    b2 = as.integer(cut(latent_b + rnorm(120, sd = 0.3), 5, labels = FALSE))
  )

  out <- engine_copula(
    data = df,
    variables = names(df),
    measurement_blocks = list(alpha = c("a1", "a2"), beta = c("b1", "b2")),
    between_block_shrink = 0.25,
    n_aug = 40,
    seed = 11
  )

  cor_used <- out$models$all$correlation
  expect_gt(abs(cor_used["a1", "a2"]), abs(cor_used["a1", "b1"]))
  expect_equal(out$models$all$between_block_shrink, 0.25)
  expect_equal(names(out$models$all$measurement_blocks), c("alpha", "beta"))
  expect_equal(out$models$all$measurement_mode, "block_copula_calibrated")
  expect_equal(out$models$all$measurement_blocks_source, "user_supplied")
  expect_equal(out$models$all$shrinkage_scope, "between_structural_families")
  expect_equal(out$models$all$score_calibration_mode, "iman_conover_family_bundle")
  expect_equal(out$models$all$score_calibration_draws, 25)
  expect_lte(
    out$models$all$structural_correlation_error_calibrated,
    out$models$all$structural_correlation_error_raw
  )
  expect_lte(
    out$models$all$structural_correlation_max_error_calibrated,
    out$models$all$structural_correlation_max_error_raw
  )
  expect_lt(out$models$all$within_family_correlation_max_change, 1e-12)
})

test_that("measurement-aware copula accepts hierarchical and covariate-only composites", {
  set.seed(14)
  n <- 100
  latent_grit <- rnorm(n)
  latent_e <- 0.3 * latent_grit + rnorm(n, sd = 0.9)
  df <- data.frame(
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

  out <- engine_copula(
    data = df,
    variables = setdiff(names(df), c("mom", "dad")),
    covariates = c("mom", "dad"),
    measurement_blocks = composite_list,
    n_aug = 40,
    seed = 15
  )

  expect_equal(names(out$models$all$measurement_blocks), c("extraversion", "grit"))
  expect_equal(sort(out$models$all$measurement_blocks$grit), sort(c("gs1", "gs2", "gs3", "gs4")))
  expect_equal(out$models$all$ignored_composites, "ses")
  expect_equal(out$models$all$measurement_mode, "block_copula_calibrated")
  expect_true(all(out$data$mom %in% df$mom))
  expect_true(all(out$data$dad %in% df$dad))
})

test_that("copula seed controls feeder noise without changing caller RNG", {
  set.seed(18)
  df <- data.frame(x = rnorm(100), y = rnorm(100))
  run <- function(ambient, seed) {
    set.seed(ambient)
    before <- .Random.seed
    out <- engine_copula(data = df, variables = c("x", "y"), n_aug = 50,
                         noise_method = "gaussian", c = 0.3, seed = seed)
    list(data = out$data, before = before, after = .Random.seed)
  }
  a <- run(101, 42)
  b <- run(202, 42)
  different <- run(101, 43)
  expect_equal(a$data, b$data)
  expect_identical(a$before, a$after)
  expect_identical(b$before, b$after)
  expect_false(identical(a$data, different$data))
})
