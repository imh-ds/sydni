test_that("discrete integer noise uses {-1, 0, +1}", {
  x <- 1:100
  y <- inject_noise(x, method = "discrete", c = 0.5)

  expect_true(all((y - x) %in% c(-1L, 0L, 1L)))
})

test_that("random noise honors user-supplied r through add_noise", {
  x <- data.frame(a = rep(10, 200))
  y <- add_noise(x, variables = "a", method = "random", r = 3, round_int = FALSE)

  expect_true(all(y$a >= 7))
  expect_true(all(y$a <= 13))
})

test_that("continuous auto noise controls remain numeric", {
  x <- rnorm(200)
  y <- inject_noise(x, method = "gaussian", c = "auto", f = 0.25, round_int = FALSE)

  expect_type(y, "double")
  expect_equal(length(y), length(x))
})
