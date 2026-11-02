# Real-data companion benchmark for the first two historical SYDNI experiments.
#
# Unlike the simulation benchmark, these datasets do not provide known truth.
# The target is therefore empirical perturbation: how much each augmentation
# method changes the observed-only analysis when synthetic rows are pooled with
# observed rows.

user_lib <- Sys.getenv("R_LIBS_USER")
local_lib <- "local_r_lib"
.libPaths(c(user_lib, local_lib, .libPaths()))

set.seed(20260612)

source_r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in source_r_files) {
  source(file, local = globalenv())
}

dir.create("results/benchmarks", recursive = TRUE, showWarnings = FALSE)

# Persistent cmdstanr model cache - see the matching comment in
# run_inferential_safety_pilot.R. Avoids recompiling on every call and avoids
# writing a fresh unsigned .exe to Temp each time (which repeatedly triggered
# Windows Defender during BRM testing).
stan_cache_dir <- normalizePath(file.path(getwd(), "stan_cache"), mustWork = FALSE)
dir.create(stan_cache_dir, recursive = TRUE, showWarnings = FALSE)
options(cmdstanr_write_stan_file_dir = stan_cache_dir)

output_stem <- Sys.getenv("SYDNI_BENCH_OUTPUT_STEM", "realdata_experiment_pilot")
method_mode <- Sys.getenv("SYDNI_BENCH_METHOD_MODE", "all")
reps <- as.integer(Sys.getenv("SYDNI_BENCH_REPS", "1"))
sydni_noise_method <- "discrete"
sydni_noise_p <- 0.80
sydni_shrinkage_k <- 2
sydni_brm_n_aug_cap <- 10L

known_method_modes <- c("all", "non_brm", "brm_only")
if (!method_mode %in% known_method_modes) {
  warning(
    "SYDNI_BENCH_METHOD_MODE = \"", method_mode, "\" is not one of ",
    paste(sprintf('"%s"', known_method_modes), collapse = ", "),
    "; falling back to the non-BRM method set. If this was unintentional, ",
    "check whether the environment variable was left over from a previous run."
  )
}

# Env vars are read once above and persist for the whole run. Printing the
# resolved values here makes it obvious when a leftover value from a prior
# invocation (e.g. SYDNI_BENCH_METHOD_MODE left at "non_brm" in the same
# shell session) is silently changing this run's method coverage - this is
# how two runs intended to differ (with vs. without BRM) previously produced
# byte-identical output files with BRM missing from both.
message(
  "SYDNI_BENCH_OUTPUT_STEM = \"", output_stem, "\"\n",
  "SYDNI_BENCH_METHOD_MODE = \"", method_mode, "\"\n",
  "SYDNI_BENCH_REPS = ", reps
)

starts_with <- function(nms, prefix) {
  nms[startsWith(nms, prefix)]
}

median_impute <- function(data) {
  data[] <- lapply(data, function(x) {
    if (is.numeric(x) || is.integer(x)) {
      replacement <- stats::median(x, na.rm = TRUE)
      if (is.na(replacement)) replacement <- 0
      x[is.na(x)] <- replacement
      x
    } else {
      tab <- sort(table(x), decreasing = TRUE)
      replacement <- if (length(tab)) names(tab)[1] else ""
      x[is.na(x) | x == ""] <- replacement
      x
    }
  })
  data
}

row_mean <- function(data, vars) {
  vars <- intersect(vars, names(data))
  if (!length(vars)) {
    return(rep(NA_real_, nrow(data)))
  }
  rowMeans(data[, vars, drop = FALSE], na.rm = TRUE)
}

add_exp1_composites <- function(data) {
  data$ses <- row_mean(data, c("mom", "dad"))
  data$material <- row_mean(data, paste0("mss_", 1:4))
  data$dissatis <- row_mean(data, paste0("dcs_", 1:4))
  data$hardwork <- row_mean(data, paste0("hssb_", 1:2))
  data$ability <- row_mean(data, paste0("hssb_", 3:4))
  data$merit <- row_mean(data, c("hardwork", "ability"))
  data$socialconnections <- row_mean(data, paste0("hssb_", 5:6))
  data$privilege <- row_mean(data, paste0("hssb_", 7:8))
  data$opportunity <- row_mean(data, c("socialconnections", "privilege"))
  data$luck <- row_mean(data, paste0("hssb_", 9:10))
  data$fate <- row_mean(data, paste0("hssb_", 11:12))
  data$chance <- row_mean(data, c("luck", "fate"))
  data
}

add_exp2_composites <- function(data) {
  data$perceive_inequality <- row_mean(data, paste0("cc1_", 1:8))
  data$accept_inequality <- row_mean(data, paste0("sei_", c(2, 3, 5)))
  data$effort <- row_mean(data, sprintf("hsbm_%02d", 1:2))
  data$ability <- row_mean(data, sprintf("hsbm_%02d", 3:4))
  data$connections <- row_mean(data, sprintf("hsbo_%02d", 1:2))
  data$privilege <- row_mean(data, sprintf("hsbo_%02d", 3:4))
  data$luck <- row_mean(data, sprintf("hsbc_%02d", 1:2))
  data$fate <- row_mean(data, sprintf("hsbc_%02d", 3:4))
  data$civic_attitude <- row_mean(data, paste0("ce1_", 1:8))
  data$civic_behavior <- row_mean(data, paste0("ce2_", 1:6))
  data$merit <- row_mean(data, c("effort", "ability"))
  data$opportunity <- row_mean(data, c("connections", "privilege"))
  data$chance <- row_mean(data, c("luck", "fate"))
  data
}

prepare_exp1 <- function() {
  data <- read.csv(
    "tests/3_Experiments/Experiment_01/1_Data/cc_data_s2_data.csv",
    check.names = FALSE
  )
  data <- data[data$attention == 1 & data$Late == 0 & data$gen != 3, ]
  keep <- c(
    "ResponseId", "age", "gen", "DAD", "MOM",
    starts_with(names(data), "mss"),
    starts_with(names(data), "dcs"),
    starts_with(names(data), "hssb")
  )
  keep <- keep[!grepl("DO$", keep)]
  data <- data[, keep, drop = FALSE]
  names(data) <- tolower(names(data))
  data$responseid <- sprintf("id_%03d", seq_len(nrow(data)))
  data <- median_impute(data)

  variables <- c(paste0("mss_", 1:4), paste0("dcs_", 1:4), paste0("hssb_", 1:12))
  covariates <- c("age", "gen", "dad", "mom")
  list(
    name = "Experiment_01_materialistic_consumer",
    data = data,
    id = "responseid",
    target = "gen",
    target_value = 2,
    variables = variables,
    covariates = covariates,
    composite_fun = add_exp1_composites,
    analysis_formula = dissatis ~ material + merit + opportunity + chance + age + dad + mom,
    structure_vars = c("material", "dissatis", "merit", "opportunity", "chance"),
    n_aug = 142
  )
}

prepare_exp2 <- function() {
  data <- read.csv(
    "tests/3_Experiments/Experiment_02/1_Data/s2_data.csv",
    check.names = FALSE
  )
  drop <- c("conservative", "female", "parent_educ", "ethn")
  drop <- c(drop, starts_with(names(data), "lsb"), starts_with(names(data), "cabc"), starts_with(names(data), "babm"))
  data <- data[, setdiff(names(data), drop), drop = FALSE]
  data <- data[data$gen != 3, ]
  data <- median_impute(data)

  variables <- names(data)[
    startsWith(names(data), "cc1") |
      startsWith(names(data), "sei") |
      startsWith(names(data), "hsb") |
      startsWith(names(data), "ce")
  ]
  covariates <- c(
    "age", "dad", "mom", "gen", "int", "ses_now", "wrkexp", "source",
    "white", "asian", "hispa", "other", "politic_1", "politic_2", "politic_3"
  )
  list(
    name = "Experiment_02_civic_engagement",
    data = data,
    id = "id",
    target = "gen",
    target_value = 2,
    variables = variables,
    covariates = intersect(covariates, names(data)),
    composite_fun = add_exp2_composites,
    analysis_formula = civic_behavior ~ civic_attitude + merit + opportunity + chance + age + dad + mom,
    structure_vars = c("civic_attitude", "civic_behavior", "merit", "opportunity", "chance"),
    n_aug = 394
  )
}

safe_run <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

take_n <- function(data, n) {
  if (nrow(data) == 0) return(data)
  data[sample(seq_len(nrow(data)), size = n, replace = nrow(data) < n), , drop = FALSE]
}

augment_smote <- function(cfg, method) {
  if (!requireNamespace("smotefamily", quietly = TRUE)) {
    stop("Package 'smotefamily' is not installed.")
  }
  data <- cfg$data
  x <- data[, setdiff(names(data), c(cfg$id, cfg$target)), drop = FALSE]
  x <- x[, vapply(x, is.numeric, logical(1)), drop = FALSE]
  target <- factor(data[[cfg$target]])
  fit <- if (method == "SMOTE") {
    smotefamily::SMOTE(x, target, K = 5, dup_size = 5)
  } else {
    smotefamily::ADAS(x, target, K = 5)
  }
  syn <- fit$syn_data
  names(syn)[names(syn) == "class"] <- cfg$target
  syn[[cfg$target]] <- suppressWarnings(as.numeric(as.character(syn[[cfg$target]])))
  syn <- syn[syn[[cfg$target]] == cfg$target_value, , drop = FALSE]
  syn <- take_n(syn, cfg$n_aug)
  missing_cols <- setdiff(names(data), names(syn))
  for (col in missing_cols) syn[[col]] <- NA
  syn <- syn[, names(data), drop = FALSE]
  syn[[cfg$id]] <- sprintf("aug_%s_%03d", tolower(method), seq_len(nrow(syn)))
  median_impute(syn)
}

augment_synthpop <- function(cfg, method) {
  if (!requireNamespace("synthpop", quietly = TRUE)) {
    stop("Package 'synthpop' is not installed.")
  }
  data <- cfg$data[, setdiff(names(cfg$data), cfg$id), drop = FALSE]
  syn_method <- rep(method, ncol(data))
  names(syn_method) <- names(data)
  fit <- synthpop::syn(
    data = data,
    method = syn_method,
    k = max(1, ceiling(cfg$n_aug / max(mean(data[[cfg$target]] == cfg$target_value), 0.01))),
    print.flag = FALSE
  )
  syn <- fit$syn
  syn <- syn[syn[[cfg$target]] == cfg$target_value, , drop = FALSE]
  syn <- take_n(syn, cfg$n_aug)
  syn[[cfg$id]] <- sprintf("aug_synthpop_%s_%03d", method, seq_len(nrow(syn)))
  median_impute(syn[, names(cfg$data), drop = FALSE])
}

augment_sydni_engine <- function(cfg, engine) {
  target_data <- cfg$data[cfg$data[[cfg$target]] == cfg$target_value, , drop = FALSE]
  target_vars <- target_data[, cfg$variables, drop = FALSE]
  n_aug_engine <- if (engine == "brm") {
    min(cfg$n_aug, sydni_brm_n_aug_cap)
  } else {
    cfg$n_aug
  }
  engine_control <- switch(
    engine,
    bart = list(num_trees = 20L, num_burn = 50L, num_draws = 50L),
    factor = list(n_factors = 2L, correlation_method = "polychoric", loading_shrink = 1),
    list()
  )
  fit <- sydni(
    data = target_vars,
    variables = cfg$variables,
    covariates = NULL,
    engine = engine,
    n_aug = n_aug_engine,
    noise_method = sydni_noise_method,
    c = "auto",
    p = sydni_noise_p,
    k = sydni_shrinkage_k,
    predictor_k = if (engine == "brm") 2L else if (engine == "bart") 5L else NULL,
    chains = if (engine == "brm") 1 else 2,
    cores = 1,
    iter = if (engine == "brm") 100 else 500,
    warmup = if (engine == "brm") 50 else 250,
    thin = 1,
    algorithm = if (engine == "brm") "meanfield" else "fullrank",
    # cmdstanr is not on CRAN, so it stays opt-in for package users (default
    # remains rstan), but for our own benchmark runs it is the only backend
    # that makes repeated BRM calls tractable - measured separately: rstan
    # plateaus around 530s/call even after warm-up, cmdstanr drops to ~11s
    # after the first compile. Ignored (harmless) for non-brm engines.
    backend = "cmdstanr",
    engine_control = engine_control,
    seed = sample.int(.Machine$integer.max, 1)
  )
  syn <- fit$data
  covariate_draws <- target_data[
    sample(seq_len(nrow(target_data)), size = nrow(syn), replace = TRUE),
    setdiff(cfg$covariates, cfg$target),
    drop = FALSE
  ]
  syn <- cbind(syn, covariate_draws)
  syn[[cfg$id]] <- sprintf("aug_sydni_%s_%03d", engine, seq_len(nrow(syn)))
  missing_cols <- setdiff(names(cfg$data), names(syn))
  for (col in missing_cols) syn[[col]] <- NA
  syn[[cfg$target]] <- cfg$target_value
  median_impute(syn[, names(cfg$data), drop = FALSE])
}

augment_method <- function(cfg, method) {
  if (method %in% c("SMOTE", "ADASYN")) {
    return(augment_smote(cfg, method))
  }
  if (method == "synthpop_sample") {
    return(augment_synthpop(cfg, "sample"))
  }
  if (method == "synthpop_cart") {
    return(augment_synthpop(cfg, "cart"))
  }
  if (method == "sydni_copula") {
    return(augment_sydni_engine(cfg, "copula"))
  }
  if (method == "sydni_bart") {
    return(augment_sydni_engine(cfg, "bart"))
  }
  if (method == "sydni_factor") {
    return(augment_sydni_engine(cfg, "factor"))
  }
  if (method == "sydni_graphical") {
    return(augment_sydni_engine(cfg, "graphical"))
  }
  if (method == "sydni_brm") {
    return(augment_sydni_engine(cfg, "brm"))
  }
  stop("Unknown method: ", method)
}

tidy_lm <- function(data, formula) {
  fit <- stats::lm(formula, data = data)
  tab <- as.data.frame(summary(fit)$coefficients)
  tab$term <- rownames(tab)
  names(tab)[1:4] <- c("estimate", "std_error", "statistic", "p_value")
  rownames(tab) <- NULL
  tab
}

cor_error <- function(observed, synthetic, vars) {
  observed <- observed[, vars, drop = FALSE]
  synthetic <- synthetic[, vars, drop = FALSE]
  observed <- observed[, vapply(observed, is.numeric, logical(1)), drop = FALSE]
  synthetic <- synthetic[, names(observed), drop = FALSE]
  if (ncol(observed) < 2 || nrow(synthetic) < 2) return(NA_real_)
  c1 <- suppressWarnings(stats::cor(observed, use = "pairwise.complete.obs"))
  c2 <- suppressWarnings(stats::cor(synthetic, use = "pairwise.complete.obs"))
  mean(abs(c1[upper.tri(c1)] - c2[upper.tri(c2)]), na.rm = TRUE)
}

near_duplicate_rate <- function(observed, synthetic, vars) {
  vars <- intersect(vars, intersect(names(observed), names(synthetic)))
  obs_key <- apply(round(observed[, vars, drop = FALSE], 6), 1, paste, collapse = "|")
  syn_key <- apply(round(synthetic[, vars, drop = FALSE], 6), 1, paste, collapse = "|")
  mean(syn_key %in% obs_key)
}

analyze_method <- function(cfg, method, rep_id) {
  observed <- cfg$data[cfg$data[[cfg$target]] == cfg$target_value, , drop = FALSE]
  observed_scored <- cfg$composite_fun(observed)
  obs_tab <- tidy_lm(observed_scored, cfg$analysis_formula)

  result <- safe_run(augment_method(cfg, method))
  failed <- inherits(result$value, "error")
  if (failed) {
    return(data.frame(
      experiment = cfg$name,
      replicate = rep_id,
      method = method,
      term = NA_character_,
      estimate_obs = NA_real_,
      estimate_pooled = NA_real_,
      se_obs = NA_real_,
      se_pooled = NA_real_,
      p_obs = NA_real_,
      p_pooled = NA_real_,
      beta_change = NA_real_,
      beta_abs_difference = NA_real_,
      effect_magnitude_change = NA_real_,
      effect_magnitude_inflation = NA_real_,
      effect_magnitude_shrinkage = NA_real_,
      se_ratio = NA_real_,
      augmentation_only_sig = NA,
      sign_flip = NA,
      cor_error = NA_real_,
      duplicate_rate = NA_real_,
      failed = TRUE,
      warning_count = length(result$warnings),
      message = conditionMessage(result$value)
    ))
  }

  synthetic <- result$value
  synthetic_scored <- cfg$composite_fun(synthetic)
  pooled_scored <- rbind(observed_scored, synthetic_scored)
  pooled_tab <- tidy_lm(pooled_scored, cfg$analysis_formula)

  merged <- merge(
    obs_tab,
    pooled_tab,
    by = "term",
    suffixes = c("_obs", "_pooled")
  )
  merged <- merged[merged$term != "(Intercept)", ]
  merged$experiment <- cfg$name
  merged$replicate <- rep_id
  merged$method <- method
  merged$beta_change <- merged$estimate_pooled - merged$estimate_obs
  merged$beta_abs_difference <- abs(merged$estimate_pooled - merged$estimate_obs)
  merged$effect_magnitude_change <- abs(merged$estimate_pooled) - abs(merged$estimate_obs)
  merged$effect_magnitude_inflation <- pmax(merged$effect_magnitude_change, 0)
  merged$effect_magnitude_shrinkage <- pmax(-merged$effect_magnitude_change, 0)
  merged$se_ratio <- merged$std_error_pooled / merged$std_error_obs
  merged$augmentation_only_sig <- merged$p_value_obs >= 0.05 & merged$p_value_pooled < 0.05
  merged$sign_flip <- sign(merged$estimate_obs) != sign(merged$estimate_pooled)
  merged$cor_error <- cor_error(observed_scored, synthetic_scored, cfg$structure_vars)
  merged$duplicate_rate <- near_duplicate_rate(observed, synthetic, cfg$variables)
  merged$failed <- FALSE
  merged$warning_count <- length(result$warnings)
  merged$message <- paste(unique(result$warnings), collapse = " | ")

  out <- merged[, c(
    "experiment", "replicate", "method", "term",
    "estimate_obs", "estimate_pooled",
    "std_error_obs", "std_error_pooled",
    "p_value_obs", "p_value_pooled",
    "beta_change", "beta_abs_difference",
    "effect_magnitude_change", "effect_magnitude_inflation",
    "effect_magnitude_shrinkage",
    "se_ratio", "augmentation_only_sig", "sign_flip",
    "cor_error", "duplicate_rate", "failed", "warning_count", "message"
  )]
  names(out)[names(out) == "std_error_obs"] <- "se_obs"
  names(out)[names(out) == "std_error_pooled"] <- "se_pooled"
  names(out)[names(out) == "p_value_obs"] <- "p_obs"
  names(out)[names(out) == "p_value_pooled"] <- "p_pooled"
  out
}

base_methods <- c(
  "SMOTE", "ADASYN", "synthpop_sample", "synthpop_cart",
  "sydni_copula", "sydni_graphical", "sydni_bart", "sydni_factor"
)

configs <- list(prepare_exp1(), prepare_exp2())

raw <- do.call(
  rbind,
  lapply(configs, function(cfg) {
    do.call(
      rbind,
      lapply(seq_len(reps), function(rep_id) {
        message("Running ", cfg$name, " replicate ", rep_id, "/", reps)
        methods <- switch(
          method_mode,
          non_brm = base_methods,
          brm_only = character(),
          all = base_methods,
          base_methods
        )
        if (
          method_mode %in% c("all", "brm_only") &&
            identical(cfg$name, "Experiment_01_materialistic_consumer")
        ) {
          methods <- c(methods, "sydni_brm")
        }
        if (!length(methods)) {
          return(NULL)
        }
        do.call(rbind, lapply(methods, function(method) {
          message("  - ", method)
          analyze_method(cfg, method, rep_id)
        }))
      })
    )
  })
)

summary <- aggregate(
  cbind(
    beta_abs_difference,
    effect_magnitude_change,
    effect_magnitude_inflation,
    effect_magnitude_shrinkage,
    se_ratio,
    augmentation_only_sig,
    sign_flip,
    cor_error,
    duplicate_rate,
    failed,
    warning_count
  ) ~ experiment + method,
  data = raw,
  FUN = function(x) mean(x, na.rm = TRUE)
)

term_summary <- aggregate(
  cbind(
    beta_abs_difference,
    effect_magnitude_change,
    effect_magnitude_inflation,
    effect_magnitude_shrinkage,
    se_ratio,
    augmentation_only_sig,
    sign_flip,
    cor_error,
    duplicate_rate,
    failed,
    warning_count
  ) ~ experiment + method + term,
  data = raw,
  FUN = function(x) mean(x, na.rm = TRUE)
)

write.csv(raw, file.path("results/benchmarks", paste0(output_stem, "_raw.csv")), row.names = FALSE)
write.csv(summary, file.path("results/benchmarks", paste0(output_stem, "_model_summary.csv")), row.names = FALSE)
write.csv(term_summary, file.path("results/benchmarks", paste0(output_stem, "_term_summary.csv")), row.names = FALSE)

print(term_summary[order(term_summary$experiment, term_summary$term, term_summary$beta_abs_difference), ], row.names = FALSE, digits = 3)
