# SYDNI inferential-safety benchmark framework
#
# This file is a simulation-first scaffold for evaluating synthetic
# augmentation methods on conservative inferential safety. It is intentionally
# not a privacy benchmark for synthpop, not a classifier-balancing benchmark for
# SMOTE/ADASYN, and not a generic synthetic-data realism benchmark for GANs.

safe_require <- function(package) {
  requireNamespace(package, quietly = TRUE)
}

null_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}

clip_round_likert <- function(x, lower = 1, upper = 5) {
  pmin(pmax(round(as.numeric(x)), lower), upper)
}

make_binary_target <- function(y, q = 0.72) {
  factor(ifelse(y > stats::quantile(y, q, na.rm = TRUE), "high", "low"))
}

# Mirrors the quantile-based Likert discretization already used for the x
# predictors in simulate_sydni_data(). Used only for the sydni_factor
# benchmark condition, where the continuous simulated y outcome must be
# discretized to satisfy engine_factor()'s validation (integer-like, factor,
# or ordered columns only). This is applied to BOTH the observed and
# augmented fit within that condition (not just as an internal synthesis
# trick), so the before/after comparison stays internally valid - but it
# means sydni_factor's row is answering a different, coarser-outcome version
# of the scenario than every other method's row. Results should be read as
# directional evidence about factor's behavior, not as directly comparable
# in magnitude to the continuous-y results for the other eight methods.
discretize_to_ordinal <- function(x, n_levels = 5) {
  as.integer(cut(
    x,
    breaks = stats::quantile(x, probs = seq(0, 1, length.out = n_levels + 1), na.rm = TRUE),
    include.lowest = TRUE,
    labels = FALSE
  ))
}

simulate_sydni_data <- function(
    scenario = c(
      "null_only",
      "sparse_effects",
      "correlated_predictors",
      "psychometric",
      "moderation_imbalance",
      "nonlinear",
      "small_cell"
    ),
    n_obs = 150,
    p = 20,
    seed = NULL
) {
  scenario <- match.arg(scenario)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  p <- max(p, 8L)
  predictor_names <- paste0("x", seq_len(p))
  group_prob <- switch(
    scenario,
    moderation_imbalance = 0.20,
    small_cell = 0.08,
    0.35
  )
  group <- factor(stats::rbinom(n_obs, 1, group_prob))
  latent <- stats::rnorm(n_obs)

  x <- switch(
    scenario,
    correlated_predictors = {
      shared <- stats::rnorm(n_obs)
      matrix(
        stats::rnorm(n_obs * p, sd = 0.75),
        nrow = n_obs,
        ncol = p
      ) + shared
    },
    psychometric = {
      f1 <- stats::rnorm(n_obs)
      f2 <- stats::rnorm(n_obs)
      loadings <- c(rep(0.85, p / 2), rep(0.75, p - p / 2))
      raw <- matrix(stats::rnorm(n_obs * p, sd = 0.75), nrow = n_obs)
      raw[, seq_len(p / 2)] <- raw[, seq_len(p / 2)] + loadings[seq_len(p / 2)] * f1
      raw[, (p / 2 + 1):p] <- raw[, (p / 2 + 1):p] + loadings[(p / 2 + 1):p] * f2
      raw
    },
    {
      matrix(stats::rnorm(n_obs * p), nrow = n_obs, ncol = p)
    }
  )

  x_likert <- apply(x, 2, function(z) {
    as.integer(cut(
      z,
      breaks = stats::quantile(z, probs = seq(0, 1, length.out = 6), na.rm = TRUE),
      include.lowest = TRUE,
      labels = FALSE
    ))
  })
  x_likert <- as.data.frame(x_likert)
  names(x_likert) <- predictor_names

  true_beta <- stats::setNames(rep(0, p), predictor_names)
  interaction_truth <- 0

  y <- switch(
    scenario,
    null_only = stats::rnorm(n_obs),
    sparse_effects = {
      true_beta[1:3] <- c(0.50, 0.35, -0.30)
      as.matrix(x_likert[, 1:3]) %*% true_beta[1:3] + stats::rnorm(n_obs, sd = 1.5)
    },
    correlated_predictors = {
      true_beta[1:3] <- c(0.45, 0.30, -0.25)
      as.matrix(x_likert[, 1:3]) %*% true_beta[1:3] + stats::rnorm(n_obs, sd = 1.5)
    },
    psychometric = {
      true_beta[1:3] <- c(0.40, 0.25, 0.20)
      as.matrix(x_likert[, 1:3]) %*% true_beta[1:3] + stats::rnorm(n_obs, sd = 1.4)
    },
    moderation_imbalance = {
      true_beta[1:2] <- c(0.45, 0.35)
      interaction_truth <- 0
      0.45 * x_likert$x1 + 0.35 * x_likert$x2 + 0.25 * as.numeric(as.character(group)) +
        interaction_truth * x_likert$x1 * as.numeric(as.character(group)) +
        stats::rnorm(n_obs, sd = 1.4)
    },
    nonlinear = {
      true_beta[1:2] <- c(0.30, 0.30)
      0.60 * (x_likert$x1 >= 4) + 0.45 * x_likert$x2 * (x_likert$x3 >= 4) +
        stats::rnorm(n_obs, sd = 1.5)
    },
    small_cell = {
      true_beta[1:2] <- c(0.45, 0.35)
      0.45 * x_likert$x1 + 0.35 * x_likert$x2 + 0.25 * as.numeric(as.character(group)) +
        stats::rnorm(n_obs, sd = 1.4)
    }
  )

  y <- as.numeric(y)
  class <- make_binary_target(y)

  data <- cbind(
    data.frame(y = y),
    x_likert,
    data.frame(group = group, class = class)
  )

  list(
    data = data,
    variables = c("y", predictor_names),
    covariates = c("group", "class"),
    truth = list(
      beta = c("(Intercept)" = NA_real_, true_beta),
      interaction_x1_group1 = interaction_truth,
      scenario = scenario
    )
  )
}

augment_bootstrap <- function(data, n_syn, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  data[sample(seq_len(nrow(data)), n_syn, replace = TRUE), , drop = FALSE]
}

augment_noise_bootstrap <- function(data, variables, n_syn, seed = NULL) {
  out <- augment_bootstrap(data, n_syn, seed = seed)
  item_vars <- intersect(variables, paste0("x", seq_len(200)))
  out[item_vars] <- lapply(out[item_vars], function(x) {
    clip_round_likert(as.numeric(x) + sample(c(-1, 0, 1), length(x), TRUE, c(0.15, 0.70, 0.15)))
  })
  out
}

augment_smotefamily <- function(data, n_syn, method = c("SMOTE", "ADASYN"), seed = NULL) {
  method <- match.arg(method)
  if (!safe_require("smotefamily")) {
    stop("Package `smotefamily` is required for SMOTE/ADASYN benchmarks.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  x <- data[, setdiff(names(data), "class"), drop = FALSE]
  x$group <- as.numeric(as.character(x$group))
  fit <- if (identical(method, "SMOTE")) {
    smotefamily::SMOTE(X = x, target = data$class, K = 5, dup_size = 3)
  } else {
    smotefamily::ADAS(X = x, target = data$class, K = 5)
  }

  out <- fit$data
  if (nrow(out) < n_syn) {
    out <- out[sample(seq_len(nrow(out)), n_syn, replace = TRUE), , drop = FALSE]
  } else {
    out <- out[sample(seq_len(nrow(out)), n_syn, replace = FALSE), , drop = FALSE]
  }

  item_vars <- grep("^x[0-9]+$", names(out), value = TRUE)
  out[item_vars] <- lapply(out[item_vars], clip_round_likert)
  out$group <- factor(round(as.numeric(out$group)), levels = levels(data$group))
  out$class <- factor(out$class, levels = levels(data$class))
  out[, names(data), drop = FALSE]
}

augment_synthpop <- function(data, n_syn, mode = c("sample", "cart"), seed = NULL) {
  mode <- match.arg(mode)
  if (!safe_require("synthpop")) {
    stop("Package `synthpop` is required for synthpop benchmarks.")
  }
  method <- if (identical(mode, "sample")) {
    c("normrank", rep("sample", ncol(data) - 1L))
  } else {
    # The first variable in the visit sequence has no predictors yet, so it
    # cannot use "cart" (rpart requires at least one predictor). Give it the
    # same safe first-variable method used by the "sample" branch above.
    c("normrank", rep("cart", ncol(data) - 1L))
  }

  # `synthpop`'s cart fitting machinery internally uses "y" as a generic
  # response-variable name when it builds each column's synthesis formula.
  # `simulate_sydni_data()` names the outcome column "y", so that column
  # collides with `synthpop`'s internal name: the resulting formula ends up
  # with zero predictor columns ("problem with term 1 in model.matrix: no
  # columns are assigned"), and CART fitting then fails downstream with
  # "result would be too long a vector". Renaming the outcome column before
  # calling `syn()`, then restoring the original name, avoids the collision.
  original_names <- names(data)
  y_pos <- which(original_names == "y")
  if (length(y_pos) == 1L) {
    names(data)[y_pos] <- ".sydni_outcome"
  }

  out <- synthpop::syn(
    data = data,
    method = method,
    k = n_syn,
    minnumlevels = if (identical(mode, "sample")) 0 else 5,
    print.flag = FALSE,
    seed = null_coalesce(seed, "sample")
  )$syn
  names(out) <- original_names
  names(data) <- original_names

  item_vars <- grep("^x[0-9]+$", names(out), value = TRUE)
  out[item_vars] <- lapply(out[item_vars], clip_round_likert)
  out$group <- factor(out$group, levels = levels(data$group))
  out$class <- factor(out$class, levels = levels(data$class))
  out[, original_names, drop = FALSE]
}

augment_sydni <- function(data, variables, covariates, n_syn, engine, seed = NULL, cores = "auto") {
  sydni(
    data = data,
    variables = variables,
    covariates = covariates,
    n_aug = n_syn,
    engine = engine,
    predictor_k = 5,
    engine_control = list(num_trees = 25, num_burn = 50, num_draws = 50),
    # `backend` is ignored by every engine except "brm", so setting it here
    # unconditionally is a no-op for copula/graphical/bart/factor. cmdstanr is
    # not on CRAN, so it stays an opt-in backend for package users - but for
    # our own benchmark runs it is the only backend that makes repeated BRM
    # calls tractable (measured: rstan plateaus around 530s/call even after
    # warm-up; cmdstanr drops to ~11-18s/call after the first compile).
    # `algorithm` is left at sydni()'s own default ("fullrank"): a controlled
    # comparison against "meanfield" failed to reproduce the Pareto-k
    # warnings seen in an earlier ad hoc run even with identical seeds
    # (likely floating-point non-associativity from parallel cores
    # interacting with ADVI's iterative optimization), so there is no
    # evidence yet that "meanfield" is actually more reliable here. Revisit
    # if the full benchmark run logs a meaningful Pareto-k warning rate.
    backend = "cmdstanr",
    # `cores = "auto"` resolves to `parallel::detectCores() - 1` inside
    # engine_brm(), i.e. every single sequential call would try to claim
    # nearly every core on the machine. With `chains = 2` (sydni()'s
    # default, unchanged here) only 2 of those are ever actually used per
    # call, so "auto" is harmless when conditions run one at a time. When
    # run_inferential_safety_benchmark() dispatches conditions across
    # multiple parallel workers, each worker must instead request a small,
    # fixed core count (matching chains = 2) so N concurrent workers don't
    # each try to claim nearly the whole machine at once.
    cores = cores,
    seed = null_coalesce(seed, NA)
  )$data
}

augment_method <- function(method, sim, n_syn, seed = NULL, cores = "auto") {
  data <- sim$data
  switch(
    method,
    observed_only = data[0, , drop = FALSE],
    bootstrap = augment_bootstrap(data, n_syn, seed),
    noise_bootstrap = augment_noise_bootstrap(data, sim$variables, n_syn, seed),
    SMOTE = augment_smotefamily(data, n_syn, "SMOTE", seed),
    ADASYN = augment_smotefamily(data, n_syn, "ADASYN", seed),
    synthpop_sample = augment_synthpop(data, n_syn, "sample", seed),
    synthpop_cart = augment_synthpop(data, n_syn, "cart", seed),
    sydni_copula = augment_sydni(data, sim$variables, sim$covariates, n_syn, "copula", seed, cores = cores),
    sydni_graphical = augment_sydni(data, sim$variables, sim$covariates, n_syn, "graphical", seed, cores = cores),
    sydni_bart = augment_sydni(data, sim$variables, sim$covariates, n_syn, "bart", seed, cores = cores),
    # `engine_factor()` rejects continuous columns outright, but this
    # benchmark's simulated `y` is continuous. `data` here has already had
    # `y` discretized to an ordinal scale by run_one_benchmark() specifically
    # for this method (see discretize_to_ordinal()), so this call receives
    # the same variable set as every other method, just with a coarser `y`.
    sydni_factor = augment_sydni(data, sim$variables, sim$covariates, n_syn, "factor", seed, cores = cores),
    sydni_brm = augment_sydni(data, sim$variables, sim$covariates, n_syn, "brm", seed, cores = cores),
    stop("Unknown benchmark method: ", method)
  )
}

analysis_formula_for <- function(data) {
  predictors <- grep("^x[0-9]+$", names(data), value = TRUE)
  predictors <- predictors[seq_len(min(20L, length(predictors)))]
  stats::as.formula(paste("y ~", paste(c(predictors, "group"), collapse = " + ")))
}

fit_analysis <- function(data, include_interaction = FALSE) {
  form <- analysis_formula_for(data)
  if (include_interaction) {
    form <- stats::update(form, . ~ . + x1:group)
  }
  stats::lm(form, data = data)
}

extract_lm_table <- function(fit) {
  co <- coef(summary(fit))
  data.frame(
    term = rownames(co),
    estimate = co[, "Estimate"],
    se = co[, "Std. Error"],
    p = co[, "Pr(>|t|)"],
    conf_low = co[, "Estimate"] - 1.96 * co[, "Std. Error"],
    conf_high = co[, "Estimate"] + 1.96 * co[, "Std. Error"],
    row.names = NULL
  )
}

inferential_metrics <- function(obs_fit, aug_fit, truth, scenario) {
  obs_tab <- extract_lm_table(obs_fit)
  aug_tab <- extract_lm_table(aug_fit)
  merged <- merge(obs_tab, aug_tab, by = "term", suffixes = c("_obs", "_aug"), all = TRUE)
  true_beta <- truth$beta
  merged$truth <- unname(true_beta[merged$term])
  merged$is_null <- is.finite(merged$truth) & abs(merged$truth) < sqrt(.Machine$double.eps)
  merged$is_nonzero <- is.finite(merged$truth) & abs(merged$truth) > sqrt(.Machine$double.eps)
  discoveries <- merged$p_aug < 0.05
  false_discoveries <- merged$is_null & discoveries

  data.frame(
    # --- Off-label "pooled data as discovery" stress test -----------------
    # SYDNI's documented workflow treats the observed-only analysis as
    # primary and the augmented analysis as a secondary, cautiously
    # interpreted check; it never recommends fitting one model on
    # observed+synthetic and reading its p-value as confirmatory. These
    # three metrics deliberately test that discouraged usage anyway (the
    # same off-label framing already used for SMOTE/ADASYN in this
    # benchmark) so we can see how bad the worst case is, not to represent
    # SYDNI's intended-use performance.
    type1_rate = mean(merged$p_aug[merged$is_null] < 0.05, na.rm = TRUE),
    false_discovery_rate = if (sum(discoveries, na.rm = TRUE) > 0) {
      sum(false_discoveries, na.rm = TRUE) / sum(discoveries, na.rm = TRUE)
    } else {
      0
    },
    # Legacy, truth-agnostic version retained for continuity with earlier
    # results. It mixes true-null terms (a real "manufactured finding") with
    # true-nonzero terms (which may just be recovered statistical power) and
    # should not be read as a safety metric on its own; see the two
    # truth-split versions below.
    augmentation_only_sig_rate = mean(merged$p_obs >= 0.05 & merged$p_aug < 0.05, na.rm = TRUE),
    # --- Intended-use secondary-check safety -------------------------------
    # SYDNI's documented decision rule is: if an effect is absent in the
    # observed data but appears only after augmentation, treat it as
    # fragile/suggestive, not confirmatory. The metric that actually tests
    # whether a *correctly behaving* analyst following that rule gets misled
    # is the rate at which this happens for terms with NO true effect - a
    # manufactured finding the documented caution is specifically meant to
    # catch. This is the primary intended-use safety metric.
    false_augmentation_only_sig_rate = mean(
      (merged$p_obs >= 0.05 & merged$p_aug < 0.05)[merged$is_null],
      na.rm = TRUE
    ),
    # The same event for true-nonzero terms is a different, non-alarmist
    # story: an underpowered real effect becoming detectable after
    # augmentation. Tracked separately so it isn't counted as danger.
    power_recovery_rate = mean(
      (merged$p_obs >= 0.05 & merged$p_aug < 0.05)[merged$is_nonzero],
      na.rm = TRUE
    ),
    mean_beta_abs_error_vs_truth = mean(abs(merged$estimate_aug - merged$truth), na.rm = TRUE),
    mean_beta_abs_difference_vs_obs = mean(abs(merged$estimate_aug - merged$estimate_obs), na.rm = TRUE),
    mean_effect_magnitude_change_vs_truth = mean(abs(merged$estimate_aug) - abs(merged$truth), na.rm = TRUE),
    mean_effect_magnitude_change_vs_obs = mean(abs(merged$estimate_aug) - abs(merged$estimate_obs), na.rm = TRUE),
    mean_effect_magnitude_inflation_vs_obs = mean(
      pmax(abs(merged$estimate_aug) - abs(merged$estimate_obs), 0),
      na.rm = TRUE
    ),
    mean_effect_magnitude_shrinkage_vs_obs = mean(
      pmax(abs(merged$estimate_obs) - abs(merged$estimate_aug), 0),
      na.rm = TRUE
    ),
    # Raw-unit differences above mix true-null terms (estimate_obs is pure
    # sampling noise near 0) with true-nonzero terms (estimate_obs has real,
    # non-degenerate magnitude), so a "-0.007" change cannot be read as small
    # or large without knowing what it's a fraction of. The two metrics below
    # give a proportional answer instead of a raw-unit one.
    #
    # For true-nonzero terms only, the ratio |beta_aug| / |beta_obs| is
    # well-behaved (estimate_obs is rarely near zero for a real effect) and
    # directly answers "did the effect shrink/inflate proportionally to its
    # own size": 1 = unchanged, <1 = shrank, >1 = inflated. Guarded against
    # near-zero denominators so an unlucky noisy replicate can't produce an
    # unstable ratio.
    mean_effect_ratio_nonzero_vs_obs = mean(
      ifelse(
        merged$is_nonzero & abs(merged$estimate_obs) > 1e-6,
        abs(merged$estimate_aug) / abs(merged$estimate_obs),
        NA_real_
      ),
      na.rm = TRUE
    ),
    # For all terms including true-null ones, estimate_obs itself is not a
    # safe denominator (it can land arbitrarily close to 0 by chance for a
    # true-null term), but se_obs is always well away from zero. Expressing
    # both estimates in units of the *original* standard error gives a
    # proportional, scenario-comparable version of
    # mean_effect_magnitude_change_vs_obs that stays well-defined for null
    # terms: positive means the effect grew relative to its original
    # precision, negative means it shrank.
    mean_std_effect_change_vs_obs = mean(
      abs(merged$estimate_aug) / merged$se_obs - abs(merged$estimate_obs) / merged$se_obs,
      na.rm = TRUE
    ),
    mean_se_ratio = mean(merged$se_aug / merged$se_obs, na.rm = TRUE),
    coverage_95 = mean(merged$conf_low_aug <= merged$truth & merged$conf_high_aug >= merged$truth, na.rm = TRUE),
    sign_reversal_rate = mean(
      sign(merged$estimate_aug[merged$is_nonzero]) != sign(merged$truth[merged$is_nonzero]),
      na.rm = TRUE
    ),
    dropped_term_rate = mean(is.na(merged$estimate_aug))
  )
}

structure_metrics <- function(observed, augmented) {
  item_vars <- intersect(grep("^x[0-9]+$", names(observed), value = TRUE), names(augmented))
  cor_obs <- stats::cor(observed[, item_vars, drop = FALSE], use = "pairwise.complete.obs")
  cor_aug <- stats::cor(augmented[, item_vars, drop = FALSE], use = "pairwise.complete.obs")
  null_edges <- abs(cor_obs) < 0.10
  diag(null_edges) <- FALSE
  # Pairs the observed data already showed a meaningful (non-null)
  # correlation for. mean_abs_cor_error below is mechanically larger for any
  # engine that deliberately shrinks correlations toward a conservative
  # target - a shrinking engine doing exactly what it's designed to do will
  # always look "worse" on that metric than one that doesn't shrink at all.
  # That number alone conflates "the correlation moved" with "the
  # correlation is now wrong". The two metrics below are restricted to
  # non-null pairs and separate those two questions: cor_sign_reversal_rate
  # asks whether a real relationship's *direction* survived (0.40 -> 0.25 is
  # fine; 0.40 -> -0.10 is a real problem), and
  # mean_cor_magnitude_change_nonnull reports shrinkage/inflation as a
  # signed, neutral-by-default quantity rather than an unqualified error.
  # Near-null pairs are deliberately excluded from both, since a sign flip
  # on noise-level correlations (e.g. 0.01 -> -0.02) isn't meaningful.
  non_null_edges <- !null_edges
  diag(non_null_edges) <- FALSE
  category_error <- mean(vapply(item_vars, function(v) {
    obs_tab <- prop.table(table(factor(observed[[v]], levels = 1:5)))
    aug_tab <- prop.table(table(factor(augmented[[v]], levels = 1:5)))
    mean(abs(obs_tab - aug_tab))
  }, numeric(1)))

  data.frame(
    mean_abs_cor_error = mean(abs(cor_aug - cor_obs), na.rm = TRUE),
    null_edge_spurious_rate = mean(abs(cor_aug[null_edges]) > 0.20, na.rm = TRUE),
    cor_sign_reversal_rate = mean(
      sign(cor_aug[non_null_edges]) != sign(cor_obs[non_null_edges]),
      na.rm = TRUE
    ),
    mean_cor_magnitude_change_nonnull = mean(
      abs(cor_aug[non_null_edges]) - abs(cor_obs[non_null_edges]),
      na.rm = TRUE
    ),
    likert_category_error = category_error,
    subgroup_prop_error = sum(abs(prop.table(table(augmented$group)) - prop.table(table(observed$group)))) / 2
  )
}

diagnostic_metrics <- function(augmented, elapsed, error = NULL, warning_count = 0) {
  data.frame(
    failed = !is.null(error),
    error_message = if (is.null(error)) NA_character_ else conditionMessage(error),
    error_time = if (is.null(error)) NA_character_ else format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    warning_count = warning_count,
    runtime_sec = elapsed,
    duplicate_rate = if (is.null(error)) mean(duplicated(augmented)) else NA_real_,
    constant_column_rate = if (is.null(error)) {
      mean(vapply(augmented, function(x) length(unique(stats::na.omit(x))) <= 1, logical(1)))
    } else {
      NA_real_
    }
  )
}

empty_inferential_metrics <- function() {
  data.frame(
    type1_rate = NA_real_,
    false_discovery_rate = NA_real_,
    augmentation_only_sig_rate = NA_real_,
    false_augmentation_only_sig_rate = NA_real_,
    power_recovery_rate = NA_real_,
    mean_beta_abs_error_vs_truth = NA_real_,
    mean_beta_abs_difference_vs_obs = NA_real_,
    mean_effect_magnitude_change_vs_truth = NA_real_,
    mean_effect_magnitude_change_vs_obs = NA_real_,
    mean_effect_magnitude_inflation_vs_obs = NA_real_,
    mean_effect_magnitude_shrinkage_vs_obs = NA_real_,
    mean_effect_ratio_nonzero_vs_obs = NA_real_,
    mean_std_effect_change_vs_obs = NA_real_,
    mean_se_ratio = NA_real_,
    coverage_95 = NA_real_,
    sign_reversal_rate = NA_real_,
    dropped_term_rate = NA_real_
  )
}

empty_structure_metrics <- function() {
  data.frame(
    mean_abs_cor_error = NA_real_,
    null_edge_spurious_rate = NA_real_,
    cor_sign_reversal_rate = NA_real_,
    mean_cor_magnitude_change_nonnull = NA_real_,
    likert_category_error = NA_real_,
    subgroup_prop_error = NA_real_
  )
}

# run_one_benchmark() is a thin backstop wrapper around
# run_one_benchmark_impl(): the inner tryCatch()es below already give
# richer, elapsed-time-aware failure rows for the specific failure points
# seen so far (data simulation is not one of them, deliberately - see the
# comment above simulate_sydni_data() below). This outer layer exists so
# that ANY future error type anywhere in the per-condition path - not just
# ones already anticipated - still becomes one recorded failed row instead
# of crashing the whole worker (and, via checkForRemoteErrors(), the whole
# parallel run). Never remove this without keeping some equivalent
# top-level guard.
run_one_benchmark <- function(
    scenario, n_obs, aug_ratio, method, replicate, p = 20, seed = 1000, cores = "auto"
) {
  tryCatch(
    run_one_benchmark_impl(scenario, n_obs, aug_ratio, method, replicate, p, seed, cores),
    error = function(e) {
      cbind(
        data.frame(
          scenario = scenario,
          n_obs = n_obs,
          aug_ratio = aug_ratio,
          method = method,
          replicate = replicate,
          analysis_target = "observed_plus_synthetic",
          outcome_discretized = identical(method, "sydni_factor")
        ),
        empty_inferential_metrics(),
        empty_structure_metrics(),
        diagnostic_metrics(NULL, elapsed = NA_real_, error = e, warning_count = 0L),
        synthetic_only_mean_abs_cor_error = NA_real_
      )
    }
  )
}

run_one_benchmark_impl <- function(
    scenario,
    n_obs,
    aug_ratio,
    method,
    replicate,
    p = 20,
    seed = 1000,
    cores = "auto"
) {
  sim_seed <- seed + replicate * 10000 + n_obs + round(100 * aug_ratio)
  sim <- simulate_sydni_data(scenario = scenario, n_obs = n_obs, p = p, seed = sim_seed)
  n_syn <- max(1L, round(n_obs * aug_ratio))
  include_interaction <- identical(scenario, "moderation_imbalance")

  # engine_factor() rejects the continuous y outcome; discretize it here so
  # that BOTH the observed and augmented fit for this method see the same
  # ordinal y (see discretize_to_ordinal()). This makes sydni_factor's row
  # answer a coarser-outcome version of the scenario than every other
  # method - flagged via outcome_discretized below, not silently absorbed.
  outcome_discretized <- identical(method, "sydni_factor")
  if (outcome_discretized) {
    sim$data$y <- discretize_to_ordinal(sim$data$y)
  }

  # error_row() standardizes the "record this condition as a graceful
  # failure" shape used below for any of: a degenerate observed-data draw
  # (e.g. small_cell can, by chance, produce a subgroup with zero members
  # at small n_obs, making `group` a constant column - lm()'s contrasts
  # then throw "contrasts can be applied only to factors with 2 or more
  # levels"), an augmentation failure, or a degenerate augmented-data fit.
  # Without this, an uncaught error here crashes the entire parallel worker
  # (and, via checkForRemoteErrors(), the whole run) instead of being
  # recorded as one failed row among many - exactly what happened for
  # sydni_factor discretization-adjacent small_cell cells during
  # confirmation-run testing.
  error_row <- function(err, elapsed, warning_count) {
    cbind(
      data.frame(
        scenario = scenario,
        n_obs = n_obs,
        aug_ratio = aug_ratio,
        method = method,
        replicate = replicate,
        analysis_target = "observed_plus_synthetic",
        outcome_discretized = outcome_discretized
      ),
      empty_inferential_metrics(),
      empty_structure_metrics(),
      diagnostic_metrics(NULL, elapsed, error = err, warning_count = warning_count),
      synthetic_only_mean_abs_cor_error = NA_real_
    )
  }

  obs_fit <- tryCatch(
    fit_analysis(sim$data, include_interaction = include_interaction),
    error = function(e) e
  )
  if (inherits(obs_fit, "error")) {
    return(error_row(obs_fit, elapsed = 0, warning_count = 0L))
  }

  warnings <- character()
  start <- proc.time()[["elapsed"]]
  result <- tryCatch(
    withCallingHandlers(
      augment_method(method, sim, n_syn, seed = sim_seed + 1L, cores = cores),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  elapsed <- proc.time()[["elapsed"]] - start

  if (inherits(result, "error")) {
    return(error_row(result, elapsed, length(warnings)))
  }

  augmented <- if (identical(method, "observed_only")) {
    sim$data
  } else {
    rbind(sim$data, result)
  }
  synthetic_only <- if (nrow(result) == 0) sim$data else result
  aug_fit <- tryCatch(
    fit_analysis(augmented, include_interaction = include_interaction),
    error = function(e) e
  )
  if (inherits(aug_fit, "error")) {
    return(error_row(aug_fit, elapsed, length(warnings)))
  }

  cbind(
    data.frame(
      scenario = scenario,
      n_obs = n_obs,
      aug_ratio = aug_ratio,
      method = method,
      replicate = replicate,
      analysis_target = "observed_plus_synthetic",
      outcome_discretized = outcome_discretized
    ),
    inferential_metrics(obs_fit, aug_fit, sim$truth, scenario),
    structure_metrics(sim$data, augmented),
    diagnostic_metrics(augmented, elapsed, warning_count = length(warnings)),
    synthetic_only_mean_abs_cor_error = structure_metrics(sim$data, synthetic_only)$mean_abs_cor_error
  )
}

run_inferential_safety_benchmark <- function(
    scenarios = c(
      "null_only",
      "sparse_effects",
      "correlated_predictors",
      "psychometric",
      "moderation_imbalance",
      "nonlinear",
      "small_cell"
    ),
    n_obs_values = c(75, 150, 300),
    aug_ratios = c(0.5, 1, 2),
    methods = c(
      "observed_only",
      "bootstrap",
      "noise_bootstrap",
      "SMOTE",
      "ADASYN",
      "synthpop_sample",
      "synthpop_cart",
      "sydni_copula",
      "sydni_graphical",
      "sydni_bart",
      "sydni_factor",
      "sydni_brm"
    ),
    reps = 50,
    p = 20,
    seed = 1000,
    # When supplied, each completed condition's row is appended to this CSV
    # immediately, and any condition already present in it (from a prior,
    # interrupted run) is skipped on start. This makes the run resumable
    # across process kills / machine shutdowns instead of losing everything
    # until the very last row finishes - the original behavior kept all
    # results in memory and only wrote output once the entire grid was done.
    checkpoint_path = NULL,
    # --- Parallel dispatch across independent grid conditions -------------
    # n_workers = 1 (default) keeps the original single-process sequential
    # path below completely unchanged, byte-for-byte identical to every
    # prior run. n_workers > 1 splits the remaining `todo` conditions across
    # a PSOCK cluster (Windows-compatible; each condition is independent, so
    # this is embarrassingly parallel). Requires checkpoint_path, since
    # results are merged back together through per-worker shard files on
    # disk rather than returned in memory from the workers.
    n_workers = 1,
    # sydni_brm dominates this benchmark's runtime (~88% of total compute
    # in the 5,040-condition pilot run) and internally uses `chains = 2`,
    # so a single call only ever uses 2 cores at once regardless of how
    # many it requests. Sequential runs can safely use cores = "auto" (see
    # augment_sydni()) because there's only ever one call in flight. Once
    # multiple workers run concurrently, each must instead request a small
    # fixed count matching chains = 2, so N workers together claim roughly
    # N * worker_cores instead of each independently requesting nearly the
    # whole machine.
    worker_cores = 2,
    # Zero-argument function each worker calls once on startup to reproduce
    # this session's .libPaths(), sourced R/ files and framework file, and
    # any options (in particular `cmdstanr_write_stan_file_dir`, the
    # persistent Stan model cache - see run_inferential_safety_pilot.R). A
    # fresh PSOCK worker starts with none of the calling session's setup,
    # so this is required whenever n_workers > 1.
    worker_setup = NULL
) {
  grid <- expand.grid(
    scenario = scenarios,
    n_obs = n_obs_values,
    aug_ratio = aug_ratios,
    method = methods,
    replicate = seq_len(reps),
    stringsAsFactors = FALSE
  )
  grid_keys <- paste(
    grid$scenario, grid$n_obs, grid$aug_ratio, grid$method, grid$replicate,
    sep = "|||"
  )

  todo <- seq_len(nrow(grid))
  if (!is.null(checkpoint_path) && file.exists(checkpoint_path)) {
    existing <- utils::read.csv(checkpoint_path, stringsAsFactors = FALSE)
    required_cols <- c("scenario", "n_obs", "aug_ratio", "method", "replicate")
    if (nrow(existing) > 0 && all(required_cols %in% names(existing))) {
      done_keys <- paste(
        existing$scenario, existing$n_obs, existing$aug_ratio,
        existing$method, existing$replicate,
        sep = "|||"
      )
      todo <- which(!(grid_keys %in% done_keys))
    }
    message(sprintf(
      "Checkpoint found at %s: %d of %d conditions already complete, %d remaining.",
      checkpoint_path, nrow(grid) - length(todo), nrow(grid), length(todo)
    ))
  }

  if (n_workers <= 1) {
    new_rows <- vector("list", length(todo))
    for (j in seq_along(todo)) {
      i <- todo[[j]]
      row <- run_one_benchmark(
        scenario = grid$scenario[[i]],
        n_obs = grid$n_obs[[i]],
        aug_ratio = grid$aug_ratio[[i]],
        method = grid$method[[i]],
        replicate = grid$replicate[[i]],
        p = p,
        seed = seed
      )
      new_rows[[j]] <- row
      if (!is.null(checkpoint_path)) {
        file_exists <- file.exists(checkpoint_path)
        utils::write.table(
          row,
          file = checkpoint_path,
          sep = ",",
          row.names = FALSE,
          col.names = !file_exists,
          append = file_exists,
          qmethod = "double"
        )
      }
    }

    return(if (!is.null(checkpoint_path)) {
      utils::read.csv(checkpoint_path, stringsAsFactors = FALSE)
    } else {
      do.call(rbind, new_rows)
    })
  }

  # --- Parallel path ------------------------------------------------------
  if (is.null(checkpoint_path)) {
    stop("`checkpoint_path` is required when `n_workers` > 1 (results are merged through per-worker shard files).")
  }
  if (is.null(worker_setup) || !is.function(worker_setup)) {
    stop("`worker_setup` must be a zero-argument function that reproduces this session's library paths, sourced files, and options on each worker.")
  }

  shard_paths <- vapply(
    seq_len(n_workers),
    function(w) paste0(checkpoint_path, ".worker", w, ".csv"),
    character(1)
  )

  # Per-worker resumability, mirroring the sequential path's checkpoint
  # behavior: if a previous parallel run was interrupted, skip conditions
  # its shard files already recorded rather than recomputing them.
  shard_done_keys <- character()
  for (sp in shard_paths) {
    if (file.exists(sp)) {
      existing_shard <- utils::read.csv(sp, stringsAsFactors = FALSE)
      if (nrow(existing_shard) > 0) {
        shard_done_keys <- c(shard_done_keys, paste(
          existing_shard$scenario, existing_shard$n_obs, existing_shard$aug_ratio,
          existing_shard$method, existing_shard$replicate,
          sep = "|||"
        ))
      }
    }
  }
  if (length(shard_done_keys) > 0) {
    still_todo <- !(grid_keys[todo] %in% shard_done_keys)
    message(sprintf(
      "Worker shards found: %d additional conditions already complete, %d remaining.",
      sum(!still_todo), sum(still_todo)
    ))
    todo <- todo[still_todo]
  }

  if (length(todo) > 0) {
    worker_id <- (seq_along(todo) - 1) %% n_workers + 1
    chunks <- split(todo, worker_id)

    cl <- parallel::makeCluster(n_workers, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, worker_setup)

    parallel::clusterApply(
      cl, seq_along(chunks),
      function(w, chunks, grid, shard_paths, p, seed, worker_cores) {
        idx <- chunks[[w]]
        shard_path <- shard_paths[[w]]
        for (i in idx) {
          row <- run_one_benchmark(
            scenario = grid$scenario[[i]],
            n_obs = grid$n_obs[[i]],
            aug_ratio = grid$aug_ratio[[i]],
            method = grid$method[[i]],
            replicate = grid$replicate[[i]],
            p = p,
            seed = seed,
            cores = worker_cores
          )
          file_exists <- file.exists(shard_path)
          utils::write.table(
            row,
            file = shard_path,
            sep = ",",
            row.names = FALSE,
            col.names = !file_exists,
            append = file_exists,
            qmethod = "double"
          )
        }
        invisible(NULL)
      },
      chunks = chunks, grid = grid, shard_paths = shard_paths,
      p = p, seed = seed, worker_cores = worker_cores
    )
  }

  # Merge existing checkpoint (if any) with every shard, dedupe by key (a
  # shard could in principle overlap with the main checkpoint if it was
  # written between the initial todo computation and now), and rewrite the
  # merged result back to checkpoint_path so every caller only ever has to
  # read one file.
  pieces <- list()
  if (file.exists(checkpoint_path)) {
    pieces[[length(pieces) + 1]] <- utils::read.csv(checkpoint_path, stringsAsFactors = FALSE)
  }
  for (sp in shard_paths) {
    if (file.exists(sp)) {
      pieces[[length(pieces) + 1]] <- utils::read.csv(sp, stringsAsFactors = FALSE)
    }
  }
  merged <- do.call(rbind, pieces)
  merged_keys <- paste(
    merged$scenario, merged$n_obs, merged$aug_ratio, merged$method, merged$replicate,
    sep = "|||"
  )
  merged <- merged[!duplicated(merged_keys), , drop = FALSE]
  utils::write.csv(merged, checkpoint_path, row.names = FALSE)
  merged
}

summarize_inferential_safety <- function(results) {
  by_cols <- c("scenario", "n_obs", "aug_ratio", "method", "outcome_discretized")
  metric_cols <- c(
    "type1_rate",
    "false_discovery_rate",
    "augmentation_only_sig_rate",
    "false_augmentation_only_sig_rate",
    "power_recovery_rate",
    "mean_beta_abs_error_vs_truth",
    "mean_beta_abs_difference_vs_obs",
    "mean_effect_magnitude_change_vs_truth",
    "mean_effect_magnitude_change_vs_obs",
    "mean_effect_magnitude_inflation_vs_obs",
    "mean_effect_magnitude_shrinkage_vs_obs",
    "mean_effect_ratio_nonzero_vs_obs",
    "mean_std_effect_change_vs_obs",
    "mean_se_ratio",
    "coverage_95",
    "sign_reversal_rate",
    "mean_abs_cor_error",
    "null_edge_spurious_rate",
    "cor_sign_reversal_rate",
    "mean_cor_magnitude_change_nonnull",
    "likert_category_error",
    "subgroup_prop_error",
    "failed",
    "warning_count",
    "runtime_sec"
  )

  groups <- split(results, results[by_cols], drop = TRUE)
  rows <- lapply(groups, function(d) {
    keys <- d[1, by_cols, drop = FALSE]
    vals <- as.data.frame(as.list(stats::setNames(
      vapply(metric_cols, function(col) {
        x <- suppressWarnings(as.numeric(d[[col]]))
        if (all(is.na(x))) {
          NA_real_
        } else {
          mean(x, na.rm = TRUE)
        }
      }, numeric(1)),
      metric_cols
    )))
    cbind(keys, vals)
  })

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}
