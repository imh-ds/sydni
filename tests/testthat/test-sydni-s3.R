make_sydni_s3_fixture <- function() {
  data <- data.frame(
    continuous = seq_len(60) / 10,
    response = rep(c(2, 3, 4, 5, 6), 12),
    group = factor(rep(c("control", "treatment"), each = 30))
  )

  sydni(
    data = data,
    variables = c("continuous", "response"),
    covariates = "group",
    engine = "copula",
    n_aug = 20,
    seed = 123
  )
}

test_that("sydni returns a self-contained classed result without changing list access", {
  out <- make_sydni_s3_fixture()

  expect_s3_class(out, "sydni")
  expect_named(out, c("models", "data", "meta", "observed_data"))
  expect_equal(nrow(out$data), 20)
  expect_equal(nrow(out$observed_data), 60)
  expect_identical(names(out$observed_data), c("continuous", "response", "group"))
  expect_s3_class(out$observed_data$group, "factor")
})

test_that("sydni print and summary expose run and reproducibility details", {
  out <- make_sydni_s3_fixture()

  expect_match(paste(capture.output(print(out)), collapse = "\n"), "Engine: copula")

  details <- summary(out)
  expect_s3_class(details, "summary.sydni")
  expect_identical(details$engine, "copula")
  expect_identical(details$n_observed, 60L)
  expect_identical(details$n_synthetic, 20L)
  expect_identical(details$seed, 123)
  expect_match(
    paste(capture.output(print(details)), collapse = "\n"),
    "Reproducibility settings"
  )
})

test_that("sydni plots observed-versus-synthetic marginal and correlation diagnostics", {
  out <- make_sydni_s3_fixture()
  plot_path <- tempfile(fileext = ".png")
  grDevices::png(plot_path)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_invisible(plot(out))
  expect_invisible(plot(out, type = "correlation"))
  expect_true(file.exists(plot_path))
})

test_that("sydni correlation plots reject requests with fewer than two eligible variables", {
  out <- make_sydni_s3_fixture()

  expect_error(
    plot(out, type = "correlation", variables = "continuous"),
    "at least two"
  )
})
