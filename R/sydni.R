#' SYDNI: synthetic data augmentation for robustness checks
#'
#' @description
#' Generates perturbed synthetic datasets for robustness and sensitivity checks.
#' SYDNI supports Bayesian regression, Gaussian copula, BART, latent-factor,
#' and sparse-graphical engines, each using a bootstrap-resampled feeder
#' dataset with controlled noise injection.
#'
#' @details
#' `sydni()` is a unified front end for multiple synthesis engines:
#'
#' - `"brm"` fits one conditional Bayesian regression model per synthesized
#'   variable by default. When `measurement_blocks` are supplied, it instead
#'   fits BRM models to structural block scores and reconstructs coherent
#'   item blocks from observed donors.
#' - `"copula"` fits subgroup-aware Gaussian copulas with empirical marginal
#'   reconstruction on noisy feeder data. With `measurement_blocks`, it
#'   restores within-family latent correlations and calibrates realized
#'   between-family score correlations toward a uniformly shrunk target.
#' - `"bart"` fits regularized Bayesian additive regression trees on noisy
#'   feeder data. With `measurement_blocks`, it preserves complete item-family
#'   bundles while calibrating between-family score correlations toward a
#'   uniformly shrunk structural target.
#' - `"factor"` fits subgroup-aware latent factor models for psychometric item
#'   batteries. With `measurement_blocks`, it preserves complete generated
#'   item-family bundles while calibrating realized between-family score
#'   correlations toward a uniformly shrunk structural target.
#' - `"graphical"` fits subgroup-aware sparse Gaussian graphical models with
#'   empirical marginal reconstruction. With `measurement_blocks`, it restores
#'   within-family latent correlations and calibrates realized between-family
#'   score correlations toward a uniformly shrunk structural target.
#'
#' SYDNI is intended for robustness, sensitivity, and workflow stress-testing.
#' It is not a general solution for recovering population truth from undersized
#' samples, and its synthetic output should be validated for each use case.
#'
#' @param data A data frame or tibble.
#' @param variables Character vector of columns in `data` to synthesize. If
#'   `NULL`, all non-covariate columns are synthesized.
#' @param covariates Character vector of columns in `data` to carry through as
#'   covariates.
#' @param engine Character. One of `"brm"`, `"copula"`, `"bart"`, `"factor"`,
#'   or `"graphical"`.
#' @param noise_method Character. One of `"discrete"`, `"truncnormal"`,
#'   `"gaussian"`, or `"random"`.
#' @param c A numeric value in `(0, 1)` or `"auto"`.
#' @param round_int Logical. If `TRUE`, integer-like synthesized values are
#'   rounded after perturbation or reconstruction.
#' @param p Numeric. Probability control for discrete noise when `c = "auto"`.
#' @param f Numeric. Variance ratio control for continuous noise when `c =
#'   "auto"`.
#' @param r Numeric. Maximum magnitude of noise when `noise_method = "random"`.
#' @param noise_bounds Optional numeric vector of length 2 used to clamp noisy
#'   feeder values.
#' @param predictor_method String. Predictor-selection rule for the BRM engine.
#' @param predictor_k Integer. Number of synthesized variables to include as
#'   predictors in each BRM conditional model.
#' @param replace Logical. Determines whether bootstrap sampling uses
#'   replacement.
#' @param n_aug Integer. Number of synthetic rows to return. When retention is
#'   requested, SYDNI may generate a larger internal candidate pool.
#' @param prop_aug Numeric. Optional proportion used to derive `n_aug`.
#' @param retain_by Optional name of a covariate used to retain one target
#'   group from a larger internally generated candidate pool.
#' @param retain_value Scalar value of `retain_by` to retain.
#' @param balance_by Optional name of a categorical covariate to balance.
#'   Currently used with `balance_target = "largest_group"`.
#' @param balance_value Scalar level of `balance_by` to augment.
#' @param balance_target Character. Currently only `"largest_group"`; SYDNI
#'   returns enough target-group rows to match the largest observed group.
#' @param families Optional named list of \pkg{brms} families for BRM outcomes.
#' @param base_sd Numeric. Base standard deviation for BRM slope priors.
#' @param sd_scale Numeric. Scale multiplier for BRM intercept priors.
#' @param k Numeric. Shrinkage exponent for BRM slope priors.
#' @param chains Numeric. Number of BRM Markov chains.
#' @param cores Integer or `"auto"`. Number of cores for BRM model fitting.
#' @param iter Numeric. BRM iterations per chain.
#' @param warmup Integer. BRM warmup iterations.
#' @param thin Integer. BRM thinning interval.
#' @param algorithm String. BRM fitting algorithm.
#' @param backend Character. One of `"rstan"` or `"cmdstanr"`.
#' @param engine_control Optional named list of engine-specific tuning
#'   parameters. Used by the `"copula"`, `"bart"`, `"factor"`, and
#'   `"graphical"` engines.
#' @param measurement_blocks Optional named list defining item blocks whose
#'   within-block correlations or measurement structure should be preserved by
#'   measurement-aware engines. Currently supported by the BRM, copula, BART,
#'   factor, and graphical engines. Also accepts a
#'   `cscore::composite_list()`-style list with `lower` and `higher` elements.
#'   Variables omitted from a supplied block list are treated as singleton
#'   families.
#' @param between_block_shrink `NULL`, `"auto"`, or a numeric value in
#'   `[0, 1]`. Automatic shrinkage is
#'   `1 / (1 + n_aug / nrow(data))`, based on the synthetic rows returned, not
#'   the size of an internal candidate pool. When
#'   `measurement_blocks` is supplied to a measurement-aware engine,
#'   between-block correlations or between-block predictor counts are shrunk by
#'   this value. If blocks are supplied and the value is omitted, automatic
#'   shrinkage is used. Supplying a shrinkage value without blocks activates
#'   all-variable shrinkage through implicit singleton families.
#' @param seed Numeric. Seed controlling internal sampling, noise injection,
#'   predictor selection, and model fitting.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{models}{Model objects used to generate the synthetic data.}
#'   \item{data}{Synthetic data generated by the chosen engine.}
#'   \item{meta}{Lightweight metadata describing the generation settings and
#'   input summary.}
#' }
#'
#' @seealso [engine_brm()], [engine_copula()], [engine_bart()],
#'   [engine_factor()], [engine_graphical()]
#'
#' @examples
#' \dontrun{
#' result <- sydni(
#'   data = grit,
#'   variables = sprintf("a%01d", seq(5)),
#'   covariates = c("gs1", "gs2"),
#'   n_aug = 100,
#'   engine = "copula",
#'   seed = 123
#' )
#'
#' head(result$data)
#' }
#'
#' @export
sydni <- function(
    data,
    variables = NULL,
    covariates = NULL,
    engine = c("brm", "copula", "bart", "factor", "graphical"),
    noise_method = c("truncnormal", "gaussian", "discrete", "random"),
    c = "auto",
    round_int = TRUE,
    p = 0.5,
    f = 0.5,
    r = 1,
    noise_bounds = NULL,
    predictor_method = c("random_information", "random", "information", "full"),
    predictor_k = NULL,
    replace = TRUE,
    n_aug = NULL,
    prop_aug = NULL,
    retain_by = NULL,
    retain_value = NULL,
    balance_by = NULL,
    balance_value = NULL,
    balance_target = "largest_group",
    families = NULL,
    base_sd = 10,
    sd_scale = 2,
    k = 2,
    chains = 2,
    cores = "auto",
    iter = 500,
    warmup = floor(iter / 2),
    thin = 2,
    algorithm = c("fullrank", "meanfield", "sampling", "pathfinder", "laplace", "fixed_param"),
    backend = c("rstan", "cmdstanr"),
    engine_control = NULL,
    measurement_blocks = NULL,
    between_block_shrink = NULL,
    seed = NA
) {
  engine <- match.arg(engine)
  noise_method <- match.arg(noise_method)
  predictor_method <- match.arg(predictor_method)
  algorithm <- match.arg(algorithm)
  backend <- match.arg(backend)
  engine_control <- validate_engine_control(engine_control)

  if (backend == "cmdstanr" && engine == "brm") {
    if (!requireNamespace("cmdstanr", quietly = TRUE)) {
      stop(
        "The 'cmdstanr' package must be installed when using backend = 'cmdstanr'."
      )
    }
    if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
      stop("CmdStan must be installed and configured for backend = 'cmdstanr'.")
    }
  }

  columns <- resolve_sydni_columns(
    data = data,
    variables = variables,
    covariates = covariates
  )
  variables <- columns$variables
  covariates <- columns$covariates

  validate_sydni_data(
    data = data,
    variables = variables,
    covariates = covariates
  )

  retention_mode <- "none"
  if (!is.null(balance_by)) {
    if (!is.null(retain_by) || !is.null(retain_value)) {
      stop("Use either `retain_by`/`retain_value` or `balance_by`/`balance_value`, not both.")
    }
    if (!is.null(n_aug) || !is.null(prop_aug)) {
      stop("Do not supply `n_aug` or `prop_aug` when `balance_by` is used; the required augmentation is calculated automatically.")
    }
    if (!identical(balance_target, "largest_group")) {
      stop("`balance_target` currently supports only \"largest_group\".")
    }
    if (length(balance_by) != 1 || !is.character(balance_by) || !balance_by %in% covariates) {
      stop("`balance_by` must name one column supplied in `covariates`.")
    }
    if (is.null(balance_value) || length(balance_value) != 1) {
      stop("`balance_value` must be a single target-group value.")
    }

    group_counts <- table(data[[balance_by]], useNA = "ifany")
    target_n_obs <- sum(values_match(data[[balance_by]], balance_value))
    if (target_n_obs == 0) {
      stop("`balance_value` was not observed in `balance_by`.")
    }
    n_aug <- as.integer(max(group_counts) - target_n_obs)
    if (n_aug < 1) {
      stop("The requested balance group is already tied for the largest observed group.")
    }
    retain_by <- balance_by
    retain_value <- balance_value
    retention_mode <- "balance_largest_group"
  } else {
    if (xor(is.null(retain_by), is.null(retain_value))) {
      stop("Supply both `retain_by` and `retain_value`, or neither.")
    }
    if (!is.null(retain_by)) {
      if (length(retain_by) != 1 || !is.character(retain_by) || !retain_by %in% covariates) {
        stop("`retain_by` must name one column supplied in `covariates`.")
      }
      if (length(retain_value) != 1) {
        stop("`retain_value` must be a single target-group value.")
      }
      retention_mode <- "target_group"
    }
  }

  if (!is.null(n_aug) && !is.null(prop_aug)) {
    stop("Specify only one of `n_aug` or `prop_aug`, not both.")
  }

  if (is.null(n_aug) && is.null(prop_aug)) {
    n_aug <- nrow(data)
  }

  if (!is.null(prop_aug)) {
    if (!is.numeric(prop_aug) || length(prop_aug) != 1 || is.na(prop_aug) || prop_aug <= 0) {
      stop("`prop_aug` must be a single numeric value greater than 0.")
    }
    n_aug <- round(nrow(data) * prop_aug)
  }

  if (!is.numeric(n_aug) || length(n_aug) != 1 || is.na(n_aug) || n_aug < 1) {
    stop("`n_aug` must be a single numeric value greater than or equal to 1.")
  }

  if (n_aug %% 1 != 0) {
    warning(
      paste0(
        "`n_aug` is not an integer; rounding to the nearest whole number: ",
        round(n_aug),
        "."
      )
    )
    n_aug <- round(n_aug)
  }

  shrinkage <- resolve_between_block_shrink(
    between_block_shrink = between_block_shrink,
    n_obs = nrow(data),
    n_retained = n_aug,
    auto_when_null = !is.null(measurement_blocks) || retention_mode != "none"
  )
  between_block_shrink_requested <- between_block_shrink
  between_block_shrink <- shrinkage$value

  effective_measurement_blocks <- measurement_blocks
  measurement_blocks_source <- if (is.null(measurement_blocks)) "none" else "user_supplied"
  if (is.null(effective_measurement_blocks) && !is.null(between_block_shrink)) {
    effective_measurement_blocks <- stats::setNames(as.list(variables), variables)
    measurement_blocks_source <- "implicit_singletons"
  }
  engine_measurement_blocks <- if (
    identical(engine, "copula") && is.null(measurement_blocks)
  ) {
    NULL
  } else {
    effective_measurement_blocks
  }

  candidate_n <- n_aug
  retention_prevalence <- NULL
  if (retention_mode != "none") {
    retention_prevalence <- mean(values_match(data[[retain_by]], retain_value))
    if (!is.finite(retention_prevalence) || retention_prevalence <= 0) {
      stop("The requested retained group was not observed in `retain_by`.")
    }
    candidate_n <- compute_retention_candidate_n(
      n_retained = n_aug,
      prevalence = retention_prevalence
    )
  }

  if (candidate_n > nrow(data) && isFALSE(replace)) {
    warning(
      "`replace` cannot be FALSE when `n_aug` exceeds the number of observed rows. Setting `replace = TRUE`."
    )
    replace <- TRUE
  }

  if (!is.null(predictor_k)) {
    if (!is.numeric(predictor_k) || length(predictor_k) != 1 || is.na(predictor_k) || predictor_k < 0) {
      stop("`predictor_k` must be a single non-negative number.")
    }
    predictor_k <- min(as.integer(round(predictor_k)), max(length(variables) - 1, 0))
  }

  data <- data |>
    dplyr::select(
      dplyr::all_of(variables),
      dplyr::all_of(covariates)
    )

  meta <- list(
    variables = variables,
    covariates = covariates,
    engine = engine,
    noise_method = noise_method,
    c = c,
    round_int = round_int,
    p = p,
    f = f,
    r = r,
    noise_bounds = noise_bounds,
    predictor_method = predictor_method,
    predictor_k = predictor_k,
    replace = replace,
    n_aug = n_aug,
    candidate_n = candidate_n,
    prop_aug = prop_aug,
    retention_mode = retention_mode,
    retain_by = retain_by,
    retain_value = retain_value,
    balance_by = balance_by,
    balance_value = balance_value,
    balance_target = if (is.null(balance_by)) NULL else balance_target,
    retention_prevalence_observed = retention_prevalence,
    effective_augmentation_ratio = n_aug / nrow(data),
    families = families,
    base_sd = base_sd,
    sd_scale = sd_scale,
    k = k,
    chains = chains,
    cores = cores,
    iter = iter,
    warmup = warmup,
    thin = thin,
    algorithm = algorithm,
    backend = backend,
    engine_control = engine_control,
    measurement_blocks = effective_measurement_blocks,
    measurement_blocks_source = measurement_blocks_source,
    between_block_shrink_requested = between_block_shrink_requested,
    between_block_shrink = between_block_shrink,
    shrinkage_source = shrinkage$source,
    shrinkage_formula = if (identical(shrinkage$source, "auto_effective_augmentation")) {
      "1 / (1 + n_retained / n_obs)"
    } else {
      NULL
    },
    seed = seed,
    input_summary = list(
      n_obs = nrow(data),
      n_variables = length(variables),
      n_covariates = length(covariates),
      column_names = names(data)
    )
  )

  synth <- switch(
    engine,
    brm = engine_brm(
      data = data,
      variables = variables,
      covariates = covariates,
      noise_method = noise_method,
      noise_bounds = noise_bounds,
      c = c,
      p = p,
      r = r,
      f = f,
      round_int = round_int,
      predictor_method = predictor_method,
      predictor_k = predictor_k,
      replace = replace,
      n_aug = candidate_n,
      effective_n_aug = n_aug,
      families = families,
      base_sd = base_sd,
      sd_scale = sd_scale,
      k = k,
      chains = chains,
      cores = cores,
      iter = iter,
      warmup = warmup,
      thin = thin,
      algorithm = algorithm,
      seed = seed,
      backend = backend,
      measurement_blocks = engine_measurement_blocks,
      between_block_shrink = between_block_shrink
    ),
    copula = engine_copula(
      data = data,
      variables = variables,
      covariates = covariates,
      noise_method = noise_method,
      noise_bounds = noise_bounds,
      c = c,
      p = p,
      r = r,
      round_int = round_int,
      n_aug = candidate_n,
      effective_n_aug = n_aug,
      engine_control = engine_control,
      measurement_blocks = engine_measurement_blocks,
      between_block_shrink = between_block_shrink,
      seed = seed
    ),
    bart = engine_bart(
      data = data,
      variables = variables,
      covariates = covariates,
      noise_method = noise_method,
      noise_bounds = noise_bounds,
      c = c,
      p = p,
      f = f,
      r = r,
      round_int = round_int,
      predictor_method = predictor_method,
      predictor_k = predictor_k,
      replace = replace,
      n_aug = candidate_n,
      effective_n_aug = n_aug,
      k = k,
      engine_control = engine_control,
      measurement_blocks = engine_measurement_blocks,
      between_block_shrink = between_block_shrink,
      seed = seed
    ),
    factor = engine_factor(
      data = data,
      variables = variables,
      covariates = covariates,
      noise_method = noise_method,
      noise_bounds = noise_bounds,
      c = c,
      p = p,
      f = f,
      r = r,
      round_int = round_int,
      n_aug = candidate_n,
      effective_n_aug = n_aug,
      k = k,
      engine_control = engine_control,
      measurement_blocks = engine_measurement_blocks,
      between_block_shrink = between_block_shrink,
      seed = seed
    ),
    graphical = engine_graphical(
      data = data,
      variables = variables,
      covariates = covariates,
      noise_method = noise_method,
      noise_bounds = noise_bounds,
      c = c,
      p = p,
      f = f,
      r = r,
      round_int = round_int,
      n_aug = candidate_n,
      effective_n_aug = n_aug,
      k = k,
      engine_control = engine_control,
      measurement_blocks = engine_measurement_blocks,
      between_block_shrink = between_block_shrink,
      seed = seed
    )
  )

  if (retention_mode != "none") {
    eligible <- which(values_match(synth$data[[retain_by]], retain_value))
    if (length(eligible) < n_aug) {
      stop(
        paste0(
          "Only ", length(eligible), " eligible synthetic rows were generated, but ",
          n_aug, " were requested. Re-run with a different seed or a more common retained group."
        )
      )
    }
    retained_indices <- with_seed(seed + 9000L, sample(eligible, n_aug, replace = FALSE))
    retained_data_raw <- synth$data[retained_indices, , drop = FALSE]

    measurement_plan <- make_measurement_score_plan(
      measurement_blocks = effective_measurement_blocks,
      variables = variables
    )
    observed_scores <- make_measurement_structural_score_frame(
      data = data,
      measurement_plan = measurement_plan
    )
    observed_cor <- estimate_brm_score_correlation(observed_scores)
    target_cor <- apply_between_block_shrinkage(
      cor_mat = observed_cor,
      between_block_shrink = between_block_shrink
    )
    retained_data <- calibrate_measurement_structural_bundles(
      data = retained_data_raw,
      measurement_plan = measurement_plan,
      target_cor = target_cor,
      seed = seed + 10000L,
      calibration_draws = 25L
    )
    calibration <- attr(retained_data, "calibration")
    attr(retained_data, "calibration") <- NULL
    synth$data <- retained_data
    meta$eligible_candidate_n <- length(eligible)
    meta$retained_n <- n_aug
    meta$post_retention_calibration <- TRUE

    retention_metadata <- list(
      retention_mode = retention_mode,
      retain_by = retain_by,
      retain_value = retain_value,
      candidate_n = candidate_n,
      eligible_candidate_n = length(eligible),
      retained_n = n_aug,
      effective_augmentation_ratio = n_aug / nrow(data),
      between_block_shrink = between_block_shrink,
      shrinkage_source = shrinkage$source,
      post_retention_calibration = TRUE,
      structural_correlation_observed = observed_cor,
      structural_correlation_target = target_cor,
      structural_correlation_generated_raw = calibration$correlation_raw,
      structural_correlation_generated = calibration$correlation,
      structural_correlation_error_raw = calibration$error_raw,
      structural_correlation_error_calibrated = calibration$error,
      structural_correlation_max_error_raw = calibration$max_error_raw,
      structural_correlation_max_error_calibrated = calibration$max_error,
      score_calibration_mode = "iman_conover_family_bundle_post_retention"
    )
    if ("all" %in% names(synth$models)) {
      synth$models$all <- utils::modifyList(synth$models$all, retention_metadata)
    } else {
      synth$models$all <- retention_metadata
    }
  }

  list(
    models = synth$models,
    data = synth$data,
    meta = meta
  )
}
