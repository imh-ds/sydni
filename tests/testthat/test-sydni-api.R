test_that("sydni validates overlapping synthesized variables and covariates", {
  df <- data.frame(a = 1:10, b = 11:20)

  expect_error(
    sydni(
      data = df,
      variables = c("a", "b"),
      covariates = "b",
      engine = "copula",
      n_aug = 10
    ),
    "both synthesized and treated as covariates"
  )
})

test_that("sydni derives n_aug from prop_aug", {
  df <- data.frame(a = rep(1:5, 10), b = rep(1:5, 10))

  out <- sydni(
    data = df,
    variables = c("a", "b"),
    engine = "copula",
    prop_aug = 0.5,
    seed = 1
  )

  expect_equal(nrow(out$data), 25)
  expect_true("meta" %in% names(out))
  expect_false("data" %in% names(out$meta))
  expect_equal(out$meta$input_summary$n_obs, 50)
})

test_that("sydni keeps n_aug and prop_aug mutually exclusive", {
  df <- data.frame(a = rep(1:5, 2), b = rep(1:5, 2))

  expect_error(
    sydni(
      data = df,
      variables = c("a", "b"),
      engine = "copula",
      n_aug = 10,
      prop_aug = 0.5
    ),
    "Specify only one"
  )
})

test_that("sydni defaults n_aug to the observed row count", {
  df <- data.frame(a = rep(1:5, 2), b = rep(1:5, 2))

  out <- sydni(
    data = df,
    variables = c("a", "b"),
    engine = "copula",
    seed = 1
  )

  expect_equal(nrow(out$data), nrow(df))
  expect_equal(out$meta$n_aug, nrow(df))
})

test_that("sydni rejects unsupported synthesized column types", {
  df <- data.frame(a = 1:5, b = I(list(1, 2, 3, 4, 5)))

  expect_error(
    sydni(
      data = df,
      variables = c("a", "b"),
      engine = "copula",
      n_aug = 5
    ),
    "Unsupported synthesized column types"
  )
})
