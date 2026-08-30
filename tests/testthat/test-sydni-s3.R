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

  printed <- paste(capture.output(print(out)), collapse = "\n")
  expect_match(printed, "Engine: copula")
  expect_match(printed, "Between-block shrinkage: none")

  details <- summary(out)
  expect_s3_class(details, "summary.sydni")
  expect_identical(details$engine, "copula")
  expect_identical(details$n_observed, 60L)
  expect_identical(details$n_synthetic, 20L)
  expect_identical(details$seed, 123)
  summary_printed <- paste(capture.output(print(details)), collapse = "\n")
  expect_match(summary_printed, "Reproducibility settings")
  expect_match(summary_printed, "Between-block shrinkage: none")
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

test_that("sydni returns editable ggplot diagnostics when requested", {
  skip_if_not_installed("ggplot2")
  out <- make_sydni_s3_fixture()

  marginal <- plot(out, backend = "ggplot", variables = "response")
  expect_s3_class(marginal, "ggplot")

  correlation <- plot(
    out,
    type = "correlation",
    backend = "ggplot",
    variables = c("continuous", "response")
  )
  expect_s3_class(correlation, "ggplot")

  expect_error(
    plot(out, backend = "ggplot"),
    "one display type"
  )
})

test_that("sydni can draw a shared marginal-plot legend", {
  plot_path <- tempfile(fileext = ".png")
  grDevices::png(plot_path)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_invisible(plot_sydni_marginal_legend())
  expect_true(file.exists(plot_path))
})

test_that("sydni lays out marginal panels in variable order", {
  expect_identical(
    sydni_marginal_layout(6L, 2L),
    rbind(c(1L, 2L), c(3L, 4L), c(5L, 6L), c(7L, 7L))
  )
})

test_that("sydni treats low-cardinality integer-like numeric variables as discrete", {
  expect_true(is_sydni_discrete_variable(c(1, 2, 3, 4, 5)))
  expect_true(is_sydni_discrete_variable(c(1, 2, 3, 4, 5) + 0))
  expect_true(is_sydni_discrete_variable(factor(c("low", "high"))))
  expect_false(is_sydni_discrete_variable(c(1.1, 2.2, 3.3, 4.4)))
  expect_false(is_sydni_discrete_variable(seq_len(11)))
})

test_that("sydni orders discrete numeric categories naturally in marginal plots", {
  expect_identical(
    sydni_discrete_levels(c(3, 1, 2), c(5, 4, 2)),
    c("1", "2", "3", "4", "5")
  )
  expect_identical(
    sydni_discrete_levels(factor(c("medium", "low"), levels = c("low", "medium", "high")), "high"),
    c("low", "medium", "high")
  )
})

test_that("sydni correlation plots reject requests with fewer than two eligible variables", {
  out <- make_sydni_s3_fixture()

  expect_error(
    plot(out, type = "correlation", variables = "continuous"),
    "at least two"
  )
})
