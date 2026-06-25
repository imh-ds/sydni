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

test_that("sydni forwards measurement-aware copula arguments", {
  df <- data.frame(
    a1 = rep(1:5, 20),
    a2 = rep(5:1, 20),
    b1 = rep(c(1:5, 5:1), 10),
    b2 = rep(c(5:1, 1:5), 10)
  )

  out <- sydni(
    data = df,
    variables = names(df),
    engine = "copula",
    measurement_blocks = list(alpha = c("a1", "a2"), beta = c("b1", "b2")),
    between_block_shrink = 0.25,
    n_aug = 20,
    seed = 3
  )

  expect_equal(nrow(out$data), 20)
  expect_equal(out$meta$between_block_shrink, 0.25)
  expect_equal(out$models$all$between_block_shrink, 0.25)
  expect_equal(names(out$models$all$measurement_blocks), c("alpha", "beta"))
})

test_that("sydni forwards measurement-aware factor arguments", {
  skip_if_not_installed("psych")

  set.seed(18)
  n <- 100
  latent_a <- rnorm(n)
  latent_b <- 0.5 * latent_a + rnorm(n, sd = 0.9)
  df <- data.frame(
    a1 = as.integer(cut(latent_a + rnorm(n, sd = 0.4), 5, labels = FALSE)),
    a2 = as.integer(cut(latent_a + rnorm(n, sd = 0.4), 5, labels = FALSE)),
    b1 = as.integer(cut(latent_b + rnorm(n, sd = 0.4), 5, labels = FALSE)),
    b2 = as.integer(cut(latent_b + rnorm(n, sd = 0.4), 5, labels = FALSE))
  )

  out <- sydni(
    data = df,
    variables = names(df),
    engine = "factor",
    measurement_blocks = list(alpha = c("a1", "a2"), beta = c("b1", "b2")),
    between_block_shrink = 0.25,
    n_aug = 20,
    seed = 3
  )

  expect_equal(nrow(out$data), 20)
  expect_equal(out$meta$between_block_shrink, 0.25)
  expect_equal(out$models$all$between_block_shrink, 0.25)
  expect_equal(out$models$all$measurement_mode, "block_one_factor")
  expect_equal(names(out$models$all$measurement_blocks), c("alpha", "beta"))
})

test_that("sydni warns about measurement arguments for unsupported engines", {
  skip_if_not_installed("dbarts")

  df <- data.frame(a = rep(1:5, 2), b = rep(5:1, 2))

  seen_warning <- FALSE
  withCallingHandlers(
    sydni(
      data = df,
      variables = names(df),
      engine = "bart",
      measurement_blocks = list(alpha = c("a", "b")),
      engine_control = list(num_trees = 5, num_burn = 5, num_draws = 5),
      n_aug = 4,
      seed = 1
    ),
    warning = function(w) {
      if (grepl("currently used only", conditionMessage(w), fixed = TRUE)) {
        seen_warning <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )
  expect_true(seen_warning)
})
