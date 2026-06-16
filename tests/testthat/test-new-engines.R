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
