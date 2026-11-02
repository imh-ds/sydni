# Study 3: empirical illustration on the package's real `grit` dataset (see
# ?grit), as opposed to Simulation Studies 1 and 2's simulated DGPs. This is
# a single-run illustration, not a replicated benchmark - there is no known
# ground truth in real data, so this framework reports structure retention
# and downstream-model stability rather than type1/power recovery against a
# known truth.
#
# Design:
# - Predictors: Extraversion, Neuroticism, Agreeableness, Openness item
#   batteries (10 items each). Conscientiousness is deliberately excluded -
#   its correlation with the Grit composite (r = 0.61 in this data) is more
#   than double any other trait's, consistent with a well-documented
#   construct-redundancy critique in the grit literature (Credeacute, Tynan
#   & Harms, 2017). This is an analyst-level judgment call about the
#   substantive model, not something SYDNI itself requires.
# - Outcome: Grit, modeled as SYDNI's hierarchical measurement_blocks
#   composite-list feature (Section 4.1) for the first time - Consistency
#   of Interest (6 items) and Perseverance of Effort (6 items) as lower-order
#   leaves, combined into one higher-order Grit structural family. Neither
#   Study 1 nor Study 2 ever exercised this feature; both used flat,
#   single-level block structures.
# - 442 complete cases (of 500) used; no imputation.
# - 17 methods, single run each (aug_ratio = 1, n_syn = n_obs): the 7
#   off-label/baseline comparators from Studies 1-2, plus each of the 5
#   SYDNI engines run as "_flat" (no measurement_blocks) and "_ma"
#   (measurement-aware), mirroring Study 2's roster exactly.
# - No natural `group` covariate exists in this data (unlike the simulated
#   studies' synthetic group split). Rather than invent a fake demographic
#   split, the comparator wrappers below only require `class` - a real,
#   substantively meaningful derived variable (median split on the Grit
#   composite, chosen over a quantile split given the ~440-observation
#   sample size, so both groups retain reasonable size) - and drop the
#   `group` dependency Studies 1-2's shared wrappers hardcode. These are
#   therefore new, Study-3-specific wrappers, not edits to the shared
#   inferential_safety_framework.R helpers (which stay untouched so as not
#   to disturb Studies 1-2's already-completed, checkpointed results).

safe_require <- function(package) {
  requireNamespace(package, quietly = TRUE)
}

null_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}

clip_round_likert <- function(x, lower = 1, upper = 5) {
  pmin(pmax(round(as.numeric(x)), lower), upper)
}

make_binary_target_median <- function(y) {
  factor(ifelse(y > stats::median(y, na.rm = TRUE), "high", "low"))
}

# --- Data preparation -------------------------------------------------------

prepare_grit_case_study_data <- function() {
  grit_env <- new.env()
  utils::data("grit", package = "sydni", envir = grit_env)
  grit <- as.data.frame(grit_env$grit)

  item_blocks <- list(
    E = paste0("e", 1:10),
    N = paste0("n", 1:10),
    A = paste0("a", 1:10),
    O = paste0("o", 1:10),
    CI = c("gs2", "gs3", "gs5", "gs7", "gs8", "gs11"),
    PE = c("gs1", "gs4", "gs6", "gs9", "gs10", "gs12")
  )
  all_items <- unlist(item_blocks, use.names = FALSE)

  complete <- grit[stats::complete.cases(grit[, all_items]), , drop = FALSE]

  grit_total <- rowMeans(complete[, c(item_blocks$CI, item_blocks$PE)])
  complete$class <- make_binary_target_median(grit_total)

  list(
    data = complete[, c(all_items, "class")],
    variables = all_items,
    covariates = "class",
    # Composite-list ("cscore::composite_list()"-style) measurement_blocks:
    # E, N, A, O remain flat leaf/final families; CI and PE combine into one
    # higher-order "Grit" family. Passed directly to sydni() for "_ma" calls.
    measurement_blocks = list(
      lower = item_blocks,
      higher = list(Grit = c("CI", "PE"))
    ),
    leaf_blocks = item_blocks,
    higher = list(Grit = c("CI", "PE"))
  )
}

# --- Construct-level scoring, mirroring make_brm_structural_score_frame's --
# --- convention exactly (mean of leaf scores for a multi-leaf family) ------

compute_leaf_scores <- function(data, leaf_blocks) {
  as.data.frame(lapply(leaf_blocks, function(items) rowMeans(data[, items, drop = FALSE], na.rm = TRUE)))
}

compute_family_scores <- function(data, leaf_blocks, higher) {
  leaf_scores <- compute_leaf_scores(data, leaf_blocks)
  referenced <- unique(unlist(higher, use.names = FALSE))
  singleton_names <- setdiff(names(leaf_blocks), referenced)
  family_scores <- leaf_scores[, singleton_names, drop = FALSE]
  for (family_name in names(higher)) {
    refs <- higher[[family_name]]
    family_scores[[family_name]] <- rowMeans(leaf_scores[, refs, drop = FALSE], na.rm = TRUE)
  }
  family_scores
}

fit_grit_construct_analysis <- function(data, leaf_blocks, higher) {
  scores <- compute_family_scores(data, leaf_blocks, higher)
  names(scores) <- paste0(names(scores), "_score")
  predictors <- setdiff(names(scores), "Grit_score")
  form <- stats::as.formula(paste("Grit_score ~", paste(predictors, collapse = " + ")))
  stats::lm(form, data = scores)
}

# --- Structure metrics: leaf-level (within-scale) and family-level --------
# --- (between-construct), no ground truth available -----------------------

cronbach_alpha <- function(item_data) {
  item_data <- as.matrix(item_data)
  k <- ncol(item_data)
  if (k < 2) {
    return(NA_real_)
  }
  item_var_sum <- sum(apply(item_data, 2, stats::var, na.rm = TRUE))
  total_var <- stats::var(rowSums(item_data), na.rm = TRUE)
  if (!is.finite(total_var) || total_var <= 0) {
    return(NA_real_)
  }
  (k / (k - 1)) * (1 - item_var_sum / total_var)
}

grit_structure_metrics <- function(observed, augmented, leaf_blocks, higher) {
  alpha_obs <- vapply(leaf_blocks, function(items) cronbach_alpha(observed[, items, drop = FALSE]), numeric(1))
  alpha_aug <- vapply(leaf_blocks, function(items) cronbach_alpha(augmented[, items, drop = FALSE]), numeric(1))

  within_scale_cor_error <- mean(vapply(leaf_blocks, function(items) {
    cor_obs <- stats::cor(observed[, items, drop = FALSE], use = "pairwise.complete.obs")
    cor_aug <- stats::cor(augmented[, items, drop = FALSE], use = "pairwise.complete.obs")
    diag(cor_obs) <- NA
    diag(cor_aug) <- NA
    mean(abs(cor_aug - cor_obs), na.rm = TRUE)
  }, numeric(1)), na.rm = TRUE)

  fam_obs <- compute_family_scores(observed, leaf_blocks, higher)
  fam_aug <- compute_family_scores(augmented, leaf_blocks, higher)
  cor_obs <- stats::cor(fam_obs, use = "pairwise.complete.obs")
  cor_aug <- stats::cor(fam_aug, use = "pairwise.complete.obs")
  off_diag <- upper.tri(cor_obs)

  data.frame(
    within_scale_alpha_abs_change = mean(abs(alpha_aug - alpha_obs), na.rm = TRUE),
    within_scale_cor_error = within_scale_cor_error,
    between_construct_cor_error = mean(abs(cor_aug[off_diag] - cor_obs[off_diag]), na.rm = TRUE),
    between_construct_sign_reversal_rate = mean(sign(cor_aug[off_diag]) != sign(cor_obs[off_diag]), na.rm = TRUE),
    between_construct_magnitude_change = mean(abs(cor_aug[off_diag]) - abs(cor_obs[off_diag]), na.rm = TRUE)
  )
}

# --- Downstream model comparison, no ground truth --------------------------
# Reports each term's before/after estimate, p-value, and whether
# significance status changed - the real-data analogue of Studies 1-2's
# type1/power_recovery metrics, which require a known truth this data does
# not have.

compare_lm_fits <- function(obs_fit, aug_fit) {
  extract <- function(fit) {
    co <- coef(summary(fit))
    data.frame(term = rownames(co), estimate = co[, "Estimate"], p = co[, "Pr(>|t|)"], row.names = NULL)
  }
  obs_tab <- extract(obs_fit)
  aug_tab <- extract(aug_fit)
  merged <- merge(obs_tab, aug_tab, by = "term", suffixes = c("_obs", "_aug"), all = TRUE)
  merged$sig_obs <- merged$p_obs < 0.05
  merged$sig_aug <- merged$p_aug < 0.05
  merged$newly_significant <- !merged$sig_obs & merged$sig_aug
  merged$lost_significance <- merged$sig_obs & !merged$sig_aug
  merged$sign_reversal <- sign(merged$estimate_obs) != sign(merged$estimate_aug)
  merged
}

# --- Comparator wrappers (Study-3-specific: no `group` dependency) --------

augment_bootstrap <- function(data, n_syn, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  data[sample(seq_len(nrow(data)), n_syn, replace = TRUE), , drop = FALSE]
}

augment_noise_bootstrap <- function(data, variables, n_syn, seed = NULL) {
  out <- augment_bootstrap(data, n_syn, seed = seed)
  out[variables] <- lapply(out[variables], function(x) {
    clip_round_likert(as.numeric(x) + sample(c(-1, 0, 1), length(x), TRUE, c(0.15, 0.70, 0.15)))
  })
  out
}

augment_smotefamily <- function(data, n_syn, method = c("SMOTE", "ADASYN"), seed = NULL) {
  method <- match.arg(method)
  if (!safe_require("smotefamily")) {
    stop("Package `smotefamily` is required for SMOTE/ADASYN benchmarks.")
  }
  if (!is.null(seed)) set.seed(seed)

  x <- data[, setdiff(names(data), "class"), drop = FALSE]
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

  item_vars <- setdiff(names(data), "class")
  out[item_vars] <- lapply(out[item_vars], clip_round_likert)
  out$class <- factor(out$class, levels = levels(data$class))
  out[, names(data), drop = FALSE]
}

augment_synthpop <- function(data, n_syn, mode = c("sample", "cart"), seed = NULL) {
  mode <- match.arg(mode)
  if (!safe_require("synthpop")) {
    stop("Package `synthpop` is required for synthpop benchmarks.")
  }
  # The first variable in the visit sequence has no predictors yet, so it
  # cannot use "cart" (rpart requires at least one predictor) - give it
  # "sample" regardless of mode, matching Studies 1-2's same fix. No
  # "y"-name collision risk exists here since no column is literally named
  # "y" (unlike Studies 1-2).
  method <- if (identical(mode, "sample")) {
    rep("sample", ncol(data))
  } else {
    c("sample", rep("cart", ncol(data) - 1L))
  }

  out <- synthpop::syn(
    data = data,
    method = method,
    k = n_syn,
    minnumlevels = if (identical(mode, "sample")) 0 else 5,
    print.flag = FALSE,
    seed = null_coalesce(seed, "sample")
  )$syn

  item_vars <- setdiff(names(data), "class")
  out[item_vars] <- lapply(out[item_vars], clip_round_likert)
  out$class <- factor(out$class, levels = levels(data$class))
  out[, names(data), drop = FALSE]
}

# --- SYDNI engine dispatch ---------------------------------------------------

augment_sydni_ma <- function(
    data, variables, covariates, n_syn, engine, seed = NULL, cores = "auto",
    measurement_blocks = NULL, between_block_shrink = NULL
) {
  sydni(
    data = data,
    variables = variables,
    covariates = covariates,
    n_aug = n_syn,
    engine = engine,
    predictor_k = 5,
    engine_control = list(num_trees = 25, num_burn = 50, num_draws = 50),
    backend = "cmdstanr",
    cores = cores,
    measurement_blocks = measurement_blocks,
    between_block_shrink = between_block_shrink,
    seed = null_coalesce(seed, NA)
  )$data
}

# --- Method dispatch (17 arms), mirroring Study 2's roster exactly --------

augment_method_grit <- function(method, sim, n_syn, seed = NULL, cores = "auto") {
  data <- sim$data

  if (grepl("_ma$", method)) {
    engine <- sub("^sydni_", "", sub("_ma$", "", method))
    return(augment_sydni_ma(
      data, sim$variables, sim$covariates, n_syn, engine, seed, cores,
      measurement_blocks = sim$measurement_blocks, between_block_shrink = 0.5
    ))
  }
  if (grepl("_flat$", method)) {
    engine <- sub("^sydni_", "", sub("_flat$", "", method))
    return(augment_sydni_ma(data, sim$variables, sim$covariates, n_syn, engine, seed, cores))
  }

  switch(
    method,
    observed_only = data[0, , drop = FALSE],
    bootstrap = augment_bootstrap(data, n_syn, seed),
    noise_bootstrap = augment_noise_bootstrap(data, sim$variables, n_syn, seed),
    SMOTE = augment_smotefamily(data, n_syn, "SMOTE", seed),
    ADASYN = augment_smotefamily(data, n_syn, "ADASYN", seed),
    synthpop_sample = augment_synthpop(data, n_syn, "sample", seed),
    synthpop_cart = augment_synthpop(data, n_syn, "cart", seed),
    stop("Unknown grit case-study method: ", method)
  )
}

# --- Single-run driver -------------------------------------------------------

run_one_grit_case_study <- function(method, sim, aug_ratio = 1, seed = 20260822, cores = "auto") {
  n_obs <- nrow(sim$data)
  n_syn <- max(1L, round(n_obs * aug_ratio))

  obs_fit <- fit_grit_construct_analysis(sim$data, sim$leaf_blocks, sim$higher)

  warnings <- character()
  start <- proc.time()[["elapsed"]]
  result <- tryCatch(
    withCallingHandlers(
      augment_method_grit(method, sim, n_syn, seed = seed, cores = cores),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  elapsed <- proc.time()[["elapsed"]] - start

  if (inherits(result, "error")) {
    return(list(
      summary = data.frame(
        method = method, n_obs = n_obs, n_syn = n_syn, failed = TRUE,
        warning_count = length(warnings), runtime_sec = elapsed,
        error_message = conditionMessage(result)
      ),
      coef_table = NULL,
      structure = NULL
    ))
  }

  augmented <- if (identical(method, "observed_only")) sim$data else rbind(sim$data, result)
  aug_fit <- fit_grit_construct_analysis(augmented, sim$leaf_blocks, sim$higher)

  coef_table <- cbind(method = method, compare_lm_fits(obs_fit, aug_fit))
  structure <- cbind(
    method = method,
    grit_structure_metrics(sim$data, augmented, sim$leaf_blocks, sim$higher)
  )

  list(
    summary = data.frame(
      method = method, n_obs = n_obs, n_syn = n_syn, failed = FALSE,
      warning_count = length(warnings), runtime_sec = elapsed, error_message = NA_character_
    ),
    coef_table = coef_table,
    structure = structure
  )
}
