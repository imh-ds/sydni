#' Print a SYDNI result
#'
#' @param x A `sydni` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.sydni <- function(x, ...) {
  details <- summary(x)

  cat("<sydni> Synthetic augmentation sensitivity result\n")
  cat("  Engine: ", details$engine, "\n", sep = "")
  cat(
    "  Rows: ", details$n_observed, " observed; ",
    details$n_synthetic, " synthetic (ratio ",
    format(details$effective_augmentation_ratio, trim = TRUE), ")\n",
    sep = ""
  )
  cat("  Variables: ", paste(details$variables, collapse = ", "), "\n", sep = "")
  if (length(details$covariates) > 0) {
    cat("  Covariates: ", paste(details$covariates, collapse = ", "), "\n", sep = "")
  }
  cat("  Retention: ", details$retention_mode, "\n", sep = "")
  cat("  Measurement blocks: ", details$measurement_blocks_source, "\n", sep = "")
  cat("  Between-block shrinkage: ", details$between_block_shrink, "\n", sep = "")
  invisible(x)
}

#' Summarize a SYDNI result
#'
#' @param object A `sydni` object.
#' @param ... Unused.
#' @return A `summary.sydni` object containing run and reproducibility details.
#' @export
summary.sydni <- function(object, ...) {
  meta <- object$meta
  observed_data <- object$observed_data
  n_observed <- if (is.data.frame(observed_data)) {
    nrow(observed_data)
  } else {
    meta$input_summary$n_obs
  }

  structure(
    list(
      engine = meta$engine,
      n_observed = as.integer(n_observed),
      n_synthetic = as.integer(nrow(object$data)),
      effective_augmentation_ratio = meta$effective_augmentation_ratio,
      variables = meta$variables,
      covariates = meta$covariates,
      retention_mode = meta$retention_mode,
      measurement_blocks_source = meta$measurement_blocks_source,
      between_block_shrink = meta$between_block_shrink,
      noise_method = meta$noise_method,
      predictor_method = meta$predictor_method,
      predictor_k = meta$predictor_k,
      seed = meta$seed,
      backend = meta$backend,
      algorithm = meta$algorithm,
      measurement_blocks = meta$measurement_blocks,
      shrinkage_source = meta$shrinkage_source
    ),
    class = "summary.sydni"
  )
}

#' @export
print.summary.sydni <- function(x, ...) {
  cat("<summary.sydni>\n")
  cat("  Engine: ", x$engine, "\n", sep = "")
  cat(
    "  Data: ", x$n_observed, " observed rows; ", x$n_synthetic,
    " synthetic rows\n",
    sep = ""
  )
  cat("  Variables: ", paste(x$variables, collapse = ", "), "\n", sep = "")
  if (length(x$covariates) > 0) {
    cat("  Covariates: ", paste(x$covariates, collapse = ", "), "\n", sep = "")
  }
  cat("  Retention mode: ", x$retention_mode, "\n", sep = "")
  cat("  Measurement blocks: ", x$measurement_blocks_source, "\n", sep = "")
  cat("  Between-block shrinkage: ", x$between_block_shrink, "\n", sep = "")
  cat("Reproducibility settings\n")
  cat("  Seed: ", x$seed, "\n", sep = "")
  cat("  Noise method: ", x$noise_method, "\n", sep = "")
  cat("  Predictor method: ", x$predictor_method, "\n", sep = "")
  if (!is.null(x$predictor_k)) {
    cat("  Predictor count: ", x$predictor_k, "\n", sep = "")
  }
  if (identical(x$engine, "brm")) {
    cat("  Backend: ", x$backend, "; algorithm: ", x$algorithm, "\n", sep = "")
  }
  invisible(x)
}

#' Plot SYDNI diagnostics
#'
#' @param x A `sydni` object.
#' @param type Diagnostic to draw: `"marginal"` compares observed and synthetic
#'   variable distributions, while `"correlation"` compares their correlation
#'   matrices.
#' @param variables Optional synthesized-variable names to include.
#' @param max_vars Maximum number of variables displayed when `variables` is
#'   omitted.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
plot.sydni <- function(x, type = c("marginal", "correlation"), variables = NULL, max_vars = 6L, ...) {
  type <- match.arg(type)
  observed_data <- x$observed_data
  synthetic_data <- x$data
  if (!is.data.frame(observed_data)) {
    stop("`plot.sydni()` requires `observed_data`; regenerate this result with the current sydni package.")
  }
  if (!is.data.frame(synthetic_data)) {
    stop("`sydni` object does not contain synthetic data to plot.")
  }

  available_variables <- x$meta$variables
  if (is.null(variables)) {
    if (!is.numeric(max_vars) || length(max_vars) != 1 || is.na(max_vars) || max_vars < 1) {
      stop("`max_vars` must be one positive number.")
    }
    variables <- utils::head(available_variables, as.integer(max_vars))
  }
  if (!is.character(variables) || length(variables) < 1 || any(!variables %in% available_variables)) {
    stop("`variables` must contain one or more synthesized variable names.")
  }
  if (any(!variables %in% names(observed_data)) || any(!variables %in% names(synthetic_data))) {
    stop("Selected variables must be present in both observed and synthetic data.")
  }

  if (identical(type, "marginal")) {
    plot_sydni_marginals(observed_data, synthetic_data, variables)
  } else {
    plot_sydni_correlations(observed_data, synthetic_data, variables)
  }
  invisible(x)
}

plot_sydni_marginals <- function(observed_data, synthetic_data, variables) {
  n_columns <- min(2L, length(variables))
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(ceiling(length(variables) / n_columns), n_columns))

  for (variable in variables) {
    observed <- observed_data[[variable]]
    synthetic <- synthetic_data[[variable]]
    if (is.factor(observed) || is.ordered(observed) || is.character(observed) ||
        (is.integer(observed) && length(unique(observed[!is.na(observed)])) <= 10L)) {
      plot_sydni_discrete_marginal(observed, synthetic, variable)
    } else {
      plot_sydni_continuous_marginal(observed, synthetic, variable)
    }
  }
}

plot_sydni_discrete_marginal <- function(observed, synthetic, variable) {
  levels <- union(as.character(observed), as.character(synthetic))
  observed_prop <- prop.table(table(factor(as.character(observed), levels = levels)))
  synthetic_prop <- prop.table(table(factor(as.character(synthetic), levels = levels)))
  graphics::barplot(
    rbind(observed = observed_prop, synthetic = synthetic_prop),
    beside = TRUE,
    col = c("#4C78A8", "#F58518"),
    main = variable,
    ylab = "Proportion",
    legend.text = TRUE,
    args.legend = list(x = "topright", bty = "n")
  )
}

plot_sydni_continuous_marginal <- function(observed, synthetic, variable) {
  observed <- observed[is.finite(observed)]
  synthetic <- synthetic[is.finite(synthetic)]
  if (length(observed) == 0 || length(synthetic) == 0) {
    stop("Selected continuous variables must contain finite observed and synthetic values.")
  }
  if (length(unique(observed)) < 2 || length(unique(synthetic)) < 2) {
    x_limits <- range(c(observed, synthetic))
    if (diff(x_limits) == 0) {
      x_limits <- x_limits + c(-0.5, 0.5)
    }
    graphics::plot(stats::ecdf(observed), xlim = x_limits, main = variable, xlab = variable,
                   ylab = "Empirical CDF", col = "#4C78A8")
    graphics::lines(stats::ecdf(synthetic), col = "#F58518")
  } else {
    observed_density <- stats::density(observed)
    synthetic_density <- stats::density(synthetic)
    x_limits <- range(observed_density$x, synthetic_density$x)
    y_limits <- range(observed_density$y, synthetic_density$y)
    graphics::plot(observed_density, xlim = x_limits, ylim = y_limits, main = variable,
                   xlab = variable, ylab = "Density", col = "#4C78A8")
    graphics::lines(synthetic_density, col = "#F58518")
  }
  graphics::legend("topright", legend = c("Observed", "Synthetic"),
                   col = c("#4C78A8", "#F58518"), lty = 1, bty = "n")
}

plot_sydni_correlations <- function(observed_data, synthetic_data, variables) {
  eligible <- variables[vapply(variables, function(variable) {
    observed <- observed_data[[variable]]
    synthetic <- synthetic_data[[variable]]
    supported <- (is.numeric(observed) || is.ordered(observed)) &&
      (is.numeric(synthetic) || is.ordered(synthetic))
    supported && length(unique(observed[!is.na(observed)])) > 1 &&
      length(unique(synthetic[!is.na(synthetic)])) > 1
  }, logical(1))]
  if (length(eligible) < 2) {
    stop("Correlation plots require at least two eligible numeric or ordered variables.")
  }

  observed_matrix <- vapply(eligible, function(variable) {
    value <- observed_data[[variable]]
    if (is.ordered(value)) as.numeric(value) else value
  }, numeric(nrow(observed_data)))
  synthetic_matrix <- vapply(eligible, function(variable) {
    value <- synthetic_data[[variable]]
    if (is.ordered(value)) as.numeric(value) else value
  }, numeric(nrow(synthetic_data)))
  colnames(observed_matrix) <- eligible
  colnames(synthetic_matrix) <- eligible

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 2))
  plot_sydni_correlation_heatmap(stats::cor(observed_matrix, use = "pairwise.complete.obs"), "Observed")
  plot_sydni_correlation_heatmap(stats::cor(synthetic_matrix, use = "pairwise.complete.obs"), "Synthetic")
}

plot_sydni_correlation_heatmap <- function(correlation, title) {
  n_variables <- ncol(correlation)
  colors <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
  graphics::image(
    x = seq_len(n_variables),
    y = seq_len(n_variables),
    z = correlation[, n_variables:1, drop = FALSE],
    zlim = c(-1, 1),
    col = colors,
    axes = FALSE,
    main = title,
    xlab = "",
    ylab = ""
  )
  graphics::axis(1, at = seq_len(n_variables), labels = colnames(correlation), las = 2)
  graphics::axis(2, at = seq_len(n_variables), labels = rev(colnames(correlation)), las = 2)
  graphics::box()
}
