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
