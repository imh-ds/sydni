# Simulation Study 2: measurement-aware (two-stage) calibration benchmark.
#
# Simulation Study 1 (inferential_safety_framework.R) tests SYDNI's engines
# on flat, ungrouped manifest variables (y ~ x1 + ... + x12) and never
# exercises `measurement_blocks`/`between_block_shrink` - SYDNI's documented
# core mechanism for multi-item scale data (protect within-scale item
# reliability while shrinking between-construct associations to limit false
# discovery at the construct level) has never been tested. This file is a
# separate simulation specifically for that mechanism: genuine latent
# constructs with their own indicator item batteries, and a paired
# with/without-measurement_blocks comparison for each SYDNI engine to
# isolate what the two-stage mechanism itself contributes, holding
# everything else (seed, data, engine, n_obs, aug_ratio) fixed.
#
# Reuses generic (non-measurement-specific) helpers from Simulation Study 1
# - augment_bootstrap(), augment_noise_bootstrap(), augment_smotefamily(),
# augment_synthpop(), make_binary_target(), clip_round_likert(),
# null_coalesce(), extract_lm_table(), inferential_metrics(),
# empty_inferential_metrics() - rather than duplicating already-tested
# code. Everything else here (simulation, SYDNI dispatch, construct-level
# metrics, the benchmark runner) is new.
#
# augment_smotefamily() and augment_synthpop() hardcode dependence on
# columns literally named "group" and "class" (SMOTE/ADASYN's balancing
# target and a carried-through covariate, respectively), so this
# simulation's data includes both, exactly mirroring Study 1's pattern,
# purely so those two functions can be reused unmodified.

source("inst/benchmarks/inferential_safety_framework.R")

# --- Simulation ----------------------------------------------------------

# 5 latent constructs: F1-F4 as predictors, Y as the outcome construct.
# F1/F2 have true nonzero structural paths to Y; F3/F4 are true null.
# Items are Likert-discretized (1-5), matching Study 1's item style and
# keeping every block compatible with engine_factor()'s ordinal
# requirement.
simulate_measurement_aware_data <- function(n_obs, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  item_names <- list(
    F1 = paste0("x", 1:3),
    F2 = paste0("x", 4:8),
    F3 = paste0("x", 9:12),
    F4 = paste0("x", 13:16),
    Y  = paste0("y", 1:4)
  )
  loadings <- list(F1 = 0.80, F2 = 0.75, F3 = 0.80, F4 = 0.75, Y = 0.80)

  F1 <- stats::rnorm(n_obs)
  F2 <- stats::rnorm(n_obs)
  F3 <- stats::rnorm(n_obs)
  F4 <- stats::rnorm(n_obs)

  # Which named construct carries the true nonzero effect is randomized per
  # replicate (via the already-seeded set.seed(seed) above), independently
  # of block/item naming order. This matters specifically because
  # augment_method_measurement_aware() uses predictor_method = "full" for
  # the measurement-aware BRM path, which - for this engine's between-
  # family predictor selection - deterministically takes the first
  # between_k candidates in list order (F1, F2, ...) rather than a random
  # subset, in order to keep the compiled Stan model structurally stable
  # across replicates (a ~75x runtime improvement over letting predictor
  # identity change randomly every replicate, which forces a fresh Stan
  # compile almost every call). If F1/F2 were *always* the true predictors
  # too, that determinism would silently hand sydni_brm_ma the exactly-
  # correct predictor set on every replicate for free - a real,
  # confirmed-in-testing bias having nothing to do with measurement-aware
  # calibration itself. Randomizing which name carries the true effect
  # decouples "which predictors get selected" (fixed, for cache stability)
  # from "which predictors are actually correct" (still properly random
  # per replicate), removing that bias while keeping the speed gain.
  true_beta_values <- stats::setNames(c(0.40, 0.30, 0, 0), sample(c("F1", "F2", "F3", "F4")))
  true_beta <- true_beta_values[c("F1", "F2", "F3", "F4")]
  names(true_beta) <- paste0(names(true_beta), "_score")

  y_latent <- true_beta[["F1_score"]] * F1 + true_beta[["F2_score"]] * F2 +
    true_beta[["F3_score"]] * F3 + true_beta[["F4_score"]] * F4 +
    stats::rnorm(n_obs, sd = 1.2)

  make_items <- function(factor_scores, names_vec, loading) {
    raw <- vapply(names_vec, function(nm) {
      loading * factor_scores + stats::rnorm(length(factor_scores), sd = sqrt(1 - loading^2))
    }, numeric(length(factor_scores)))
    colnames(raw) <- names_vec
    as.data.frame(apply(raw, 2, function(z) {
      as.integer(cut(
        z,
        breaks = stats::quantile(z, probs = seq(0, 1, length.out = 6), na.rm = TRUE),
        include.lowest = TRUE,
        labels = FALSE
      ))
    }))
  }

  items <- cbind(
    make_items(F1, item_names$F1, loadings$F1),
    make_items(F2, item_names$F2, loadings$F2),
    make_items(F3, item_names$F3, loadings$F3),
    make_items(F4, item_names$F4, loadings$F4),
    make_items(y_latent, item_names$Y, loadings$Y)
  )

  group <- factor(stats::rbinom(n_obs, 1, 0.35))
  class <- make_binary_target(y_latent)

  data <- cbind(items, data.frame(group = group, class = class))

  list(
    data = data,
    variables = unlist(item_names, use.names = FALSE),
    covariates = c("group", "class"),
    measurement_blocks = item_names,
    truth = list(beta = c("(Intercept)" = NA_real_, true_beta))
  )
}

# --- Construct-level analysis ---------------------------------------------

compute_construct_scores <- function(data, measurement_blocks) {
  scores <- lapply(measurement_blocks, function(items) rowMeans(data[, items, drop = FALSE]))
  as.data.frame(scores)
}

fit_construct_analysis <- function(data, measurement_blocks) {
  scores <- compute_construct_scores(data, measurement_blocks)
  names(scores) <- paste0(names(scores), "_score")
  predictors <- setdiff(names(scores), "Y_score")
  form <- stats::as.formula(paste("Y_score ~", paste(predictors, collapse = " + ")))
  stats::lm(form, data = scores)
}

# --- Measurement-aware structure metrics -----------------------------------

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

# Tests the actual claim behind measurement-aware calibration: within-scale
# item reliability/correlation should be preserved, while between-construct
# associations shrink deliberately (and should shrink without reversing
# sign or distorting arbitrarily).
measurement_structure_metrics <- function(observed, augmented, measurement_blocks) {
  alpha_obs <- vapply(
    measurement_blocks,
    function(items) cronbach_alpha(observed[, items, drop = FALSE]),
    numeric(1)
  )
  alpha_aug <- vapply(
    measurement_blocks,
    function(items) cronbach_alpha(augmented[, items, drop = FALSE]),
    numeric(1)
  )

  within_scale_cor_error <- mean(vapply(measurement_blocks, function(items) {
    if (length(items) < 2) {
      return(NA_real_)
    }
    cor_obs <- stats::cor(observed[, items, drop = FALSE], use = "pairwise.complete.obs")
    cor_aug <- stats::cor(augmented[, items, drop = FALSE], use = "pairwise.complete.obs")
    diag(cor_obs) <- NA
    diag(cor_aug) <- NA
    mean(abs(cor_aug - cor_obs), na.rm = TRUE)
  }, numeric(1)), na.rm = TRUE)

  scores_obs <- compute_construct_scores(observed, measurement_blocks)
  scores_aug <- compute_construct_scores(augmented, measurement_blocks)
  cor_obs <- stats::cor(scores_obs, use = "pairwise.complete.obs")
  cor_aug <- stats::cor(scores_aug, use = "pairwise.complete.obs")
  off_diag <- upper.tri(cor_obs)

  data.frame(
    within_scale_alpha_abs_change = mean(abs(alpha_aug - alpha_obs), na.rm = TRUE),
    within_scale_cor_error = within_scale_cor_error,
    between_construct_cor_error = mean(abs(cor_aug[off_diag] - cor_obs[off_diag]), na.rm = TRUE),
    between_construct_sign_reversal_rate = mean(
      sign(cor_aug[off_diag]) != sign(cor_obs[off_diag]), na.rm = TRUE
    ),
    between_construct_magnitude_change = mean(
      abs(cor_aug[off_diag]) - abs(cor_obs[off_diag]), na.rm = TRUE
    )
  )
}

empty_measurement_structure_metrics <- function() {
  data.frame(
    within_scale_alpha_abs_change = NA_real_,
    within_scale_cor_error = NA_real_,
    between_construct_cor_error = NA_real_,
    between_construct_sign_reversal_rate = NA_real_,
    between_construct_magnitude_change = NA_real_
  )
}

# --- Method dispatch (17 arms) ---------------------------------------------
#
# observed_only, bootstrap, noise_bootstrap, SMOTE, ADASYN, synthpop_sample,
# synthpop_cart (7 comparators, same roster as Simulation Study 1), plus
# each of the five measurement-aware-capable SYDNI engines run twice:
# "<engine>_flat" (no measurement_blocks - identical call shape to how
# every engine was called in Study 1) and "<engine>_ma" (measurement-aware,
# with measurement_blocks + between_block_shrink supplied). The flat/ma
# pair for a given engine, run on identical simulated data with identical
# seeds, isolates what measurement_blocks itself contributes.

augment_sydni_ma <- function(
    data, variables, covariates, n_syn, engine, seed = NULL, cores = "auto",
    measurement_blocks = NULL, between_block_shrink = NULL,
    predictor_method = "random_information"
) {
  sydni(
    data = data,
    variables = variables,
    covariates = covariates,
    n_aug = n_syn,
    engine = engine,
    predictor_k = 5,
    predictor_method = predictor_method,
    engine_control = list(num_trees = 25, num_burn = 50, num_draws = 50),
    backend = "cmdstanr",
    cores = cores,
    measurement_blocks = measurement_blocks,
    between_block_shrink = between_block_shrink,
    seed = null_coalesce(seed, NA)
  )$data
}

augment_method_measurement_aware <- function(method, sim, n_syn, seed = NULL, cores = "auto") {
  data <- sim$data

  if (grepl("_ma$", method)) {
    engine <- sub("^sydni_", "", sub("_ma$", "", method))
    # predictor_method is left at the package default ("random_information",
    # same as flat mode) rather than forced to "full". "full" was tried as
    # a fix for sydni_brm_ma's slow, unstable per-condition runtime, but
    # root-caused (via a direct diff of two consecutively-compiled .stan
    # files) to something unrelated to predictor selection entirely: BRM's
    # intercept prior mean is estimated from each replicate's own data and
    # baked into the generated Stan program as a literal numeric constant
    # (e.g. `normal_lpdf(Intercept | -5.03e-19, 2)`), so almost every
    # replicate compiles a textually distinct program regardless of which
    # predictors were selected. "full" therefore bought no real speed
    # benefit while introducing a genuine bias risk (deterministic
    # first-in-list-order predictor selection could align suspiciously
    # with which construct is actually true) - not worth keeping. See
    # simulate_measurement_aware_data()'s randomized true-effect
    # assignment for the (still-worthwhile, kept) fix to that bias risk
    # independent of this reversion. sydni_brm_ma's real per-condition
    # cost is close to sydni_brm_flat's, not the ~5s a stable cache would
    # have given - budget accordingly (see run_measurement_aware_pilot.R).
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
    stop("Unknown measurement-aware benchmark method: ", method)
  )
}

# --- Per-condition runner ---------------------------------------------------

# Backstop wrapper - see the identical pattern/rationale on
# run_one_benchmark() in inferential_safety_framework.R: any error type not
# already anticipated by the inner tryCatch()es below still becomes one
# recorded failed row instead of crashing the whole worker/run.
run_one_measurement_aware_benchmark <- function(
    n_obs, aug_ratio, method, replicate, seed = 1000, cores = "auto"
) {
  tryCatch(
    run_one_measurement_aware_benchmark_impl(n_obs, aug_ratio, method, replicate, seed, cores),
    error = function(e) {
      cbind(
        data.frame(n_obs = n_obs, aug_ratio = aug_ratio, method = method, replicate = replicate),
        empty_inferential_metrics(),
        empty_measurement_structure_metrics(),
        diagnostic_metrics(NULL, elapsed = NA_real_, error = e, warning_count = 0L)
      )
    }
  )
}

run_one_measurement_aware_benchmark_impl <- function(
    n_obs, aug_ratio, method, replicate, seed = 1000, cores = "auto"
) {
  sim_seed <- seed + replicate * 10000 + n_obs + round(100 * aug_ratio)
  sim <- simulate_measurement_aware_data(n_obs = n_obs, seed = sim_seed)
  n_syn <- max(1L, round(n_obs * aug_ratio))

  # error_row() standardizes the "record this condition as a graceful
  # failure" shape - see the identical comment/pattern in
  # run_one_benchmark() (inferential_safety_framework.R) for why obs_fit
  # and aug_fit, not just augmentation, need this: a degenerate simulated
  # draw can make fit_construct_analysis() throw (e.g. a constant grouping
  # column), and without a tryCatch here that crashes the entire parallel
  # worker rather than being recorded as one failed row.
  error_row <- function(err, elapsed, warning_count) {
    cbind(
      data.frame(
        n_obs = n_obs, aug_ratio = aug_ratio, method = method, replicate = replicate
      ),
      empty_inferential_metrics(),
      empty_measurement_structure_metrics(),
      diagnostic_metrics(NULL, elapsed, error = err, warning_count = warning_count)
    )
  }

  obs_fit <- tryCatch(
    fit_construct_analysis(sim$data, sim$measurement_blocks),
    error = function(e) e
  )
  if (inherits(obs_fit, "error")) {
    return(error_row(obs_fit, elapsed = 0, warning_count = 0L))
  }

  warnings <- character()
  start <- proc.time()[["elapsed"]]
  result <- tryCatch(
    withCallingHandlers(
      augment_method_measurement_aware(method, sim, n_syn, seed = sim_seed + 1L, cores = cores),
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

  augmented <- if (identical(method, "observed_only")) sim$data else rbind(sim$data, result)
  aug_fit <- tryCatch(
    fit_construct_analysis(augmented, sim$measurement_blocks),
    error = function(e) e
  )
  if (inherits(aug_fit, "error")) {
    return(error_row(aug_fit, elapsed, length(warnings)))
  }

  cbind(
    data.frame(
      n_obs = n_obs, aug_ratio = aug_ratio, method = method, replicate = replicate
    ),
    inferential_metrics(obs_fit, aug_fit, sim$truth, "measurement_aware"),
    measurement_structure_metrics(sim$data, augmented, sim$measurement_blocks),
    diagnostic_metrics(augmented, elapsed, warning_count = length(warnings))
  )
}

# --- Grid runner (checkpointed, optionally parallel) ------------------------
#
# Mirrors run_inferential_safety_benchmark()'s checkpoint/resume and
# n_workers/worker_cores/worker_setup parallel dispatch exactly - see that
# function's comments in inferential_safety_framework.R for the full
# rationale. Kept as a separate copy (not shared code) since the two
# studies' per-condition runners take different arguments (no `scenario`
# or `p` dimension here).
run_measurement_aware_benchmark <- function(
    n_obs_values = c(75, 150, 300),
    aug_ratios = c(0.5, 1, 2),
    methods,
    reps = 10,
    seed = 1000,
    checkpoint_path = NULL,
    n_workers = 1,
    worker_cores = 2,
    worker_setup = NULL
) {
  grid <- expand.grid(
    n_obs = n_obs_values,
    aug_ratio = aug_ratios,
    method = methods,
    replicate = seq_len(reps),
    stringsAsFactors = FALSE
  )
  grid_keys <- paste(grid$n_obs, grid$aug_ratio, grid$method, grid$replicate, sep = "|||")

  todo <- seq_len(nrow(grid))
  if (!is.null(checkpoint_path) && file.exists(checkpoint_path)) {
    existing <- utils::read.csv(checkpoint_path, stringsAsFactors = FALSE)
    required_cols <- c("n_obs", "aug_ratio", "method", "replicate")
    if (nrow(existing) > 0 && all(required_cols %in% names(existing))) {
      done_keys <- paste(existing$n_obs, existing$aug_ratio, existing$method, existing$replicate, sep = "|||")
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
      row <- run_one_measurement_aware_benchmark(
        n_obs = grid$n_obs[[i]], aug_ratio = grid$aug_ratio[[i]],
        method = grid$method[[i]], replicate = grid$replicate[[i]], seed = seed
      )
      new_rows[[j]] <- row
      if (!is.null(checkpoint_path)) {
        file_exists <- file.exists(checkpoint_path)
        utils::write.table(
          row, file = checkpoint_path, sep = ",", row.names = FALSE,
          col.names = !file_exists, append = file_exists, qmethod = "double"
        )
      }
    }
    return(if (!is.null(checkpoint_path)) {
      utils::read.csv(checkpoint_path, stringsAsFactors = FALSE)
    } else {
      do.call(rbind, new_rows)
    })
  }

  if (is.null(checkpoint_path)) {
    stop("`checkpoint_path` is required when `n_workers` > 1 (results are merged through per-worker shard files).")
  }
  if (is.null(worker_setup) || !is.function(worker_setup)) {
    stop("`worker_setup` must be a zero-argument function that reproduces this session's library paths, sourced files, and options on each worker.")
  }

  shard_paths <- vapply(seq_len(n_workers), function(w) paste0(checkpoint_path, ".worker", w, ".csv"), character(1))

  shard_done_keys <- character()
  for (sp in shard_paths) {
    if (file.exists(sp)) {
      existing_shard <- utils::read.csv(sp, stringsAsFactors = FALSE)
      if (nrow(existing_shard) > 0) {
        shard_done_keys <- c(shard_done_keys, paste(
          existing_shard$n_obs, existing_shard$aug_ratio, existing_shard$method, existing_shard$replicate,
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
      function(w, chunks, grid, shard_paths, seed, worker_cores) {
        idx <- chunks[[w]]
        shard_path <- shard_paths[[w]]
        for (i in idx) {
          row <- run_one_measurement_aware_benchmark(
            n_obs = grid$n_obs[[i]], aug_ratio = grid$aug_ratio[[i]],
            method = grid$method[[i]], replicate = grid$replicate[[i]],
            seed = seed, cores = worker_cores
          )
          file_exists <- file.exists(shard_path)
          utils::write.table(
            row, file = shard_path, sep = ",", row.names = FALSE,
            col.names = !file_exists, append = file_exists, qmethod = "double"
          )
        }
        invisible(NULL)
      },
      chunks = chunks, grid = grid, shard_paths = shard_paths, seed = seed, worker_cores = worker_cores
    )
  }

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
  merged_keys <- paste(merged$n_obs, merged$aug_ratio, merged$method, merged$replicate, sep = "|||")
  merged <- merged[!duplicated(merged_keys), , drop = FALSE]
  utils::write.csv(merged, checkpoint_path, row.names = FALSE)
  merged
}
