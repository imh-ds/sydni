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
