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
})
