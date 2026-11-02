# Run a pilot version of the SYDNI inferential-safety benchmark.
#
# This pilot is intentionally smaller than the manuscript-scale benchmark.
# It is meant to validate the benchmark machinery, expose failure modes, and
# provide preliminary summaries before running hundreds of replications.

user_lib <- Sys.getenv("R_LIBS_USER")
local_lib <- "local_r_lib"
.libPaths(c(user_lib, local_lib, .libPaths()))

invisible(lapply(list.files("R", full.names = TRUE), source))
source("inst/benchmarks/inferential_safety_framework.R")

dir.create("results/benchmarks", recursive = TRUE, showWarnings = FALSE)

# Persistent cmdstanr model cache: without this, cmdstanr writes each
# compiled model to a fresh tempdir() every run, which (a) forces a full
# recompile of every submodel on every call regardless of whether the same
# model structure was already compiled earlier, and (b) creates a brand-new
# unsigned .exe in Temp on every compile, which is exactly the pattern that
# triggered repeated Windows Defender blocks during BRM testing. Pointing
# cmdstanr at a fixed, reused directory (excluded from AV scanning) fixes
# both: identical model structures reuse the already-compiled, already-
# scanned binary instead of writing a new one each time.
stan_cache_dir <- normalizePath(file.path(getwd(), "stan_cache"), mustWork = FALSE)
dir.create(stan_cache_dir, recursive = TRUE, showWarnings = FALSE)
options(cmdstanr_write_stan_file_dir = stan_cache_dir)

# Parallel dispatch: sydni_brm accounts for the large majority of this
# benchmark's runtime and internally caps itself to 2 cores per call
# (chains = 2), so running multiple independent grid conditions
# concurrently is close to embarrassingly parallel rather than fighting
# over the same cores. n_workers = 1 falls back to the original single-
# process sequential path untouched. worker_setup() reproduces this
# script's own startup (library paths, sourced files, the persistent Stan
# model cache option) on each fresh PSOCK worker, which otherwise starts
# with none of it.
n_workers <- 4 # leaves 4 of 12 physical cores free rather than claiming all of them
worker_cores <- 2
# A fresh PSOCK worker does not necessarily start in this project's working
# directory, so every path below is either absolute or set relative to an
# explicit setwd() first - relying on relative paths alone (as the rest of
# this script can, since it always runs from the project root) would be
# fragile here.
# worker_setup is built with local() rather than referencing a plain global
# variable: parallel's PSOCK workers do not receive the master's global
# environment when a function is shipped to them via clusterCall(), only a
# same-process reference that resolves to an empty environment on the
# worker side. A free variable like `project_root` would fail there with
# "object not found". Binding it via local() instead puts it in the
# function's own closure environment, which parallel DOES serialize fully.
worker_setup <- local({
  project_root <- getwd()
  function() {
    setwd(project_root)
    .libPaths(c(
      Sys.getenv("R_LIBS_USER"),
      "local_r_lib",
      .libPaths()
    ))
    invisible(lapply(list.files("R", full.names = TRUE), source))
    source("inst/benchmarks/inferential_safety_framework.R")
    options(cmdstanr_write_stan_file_dir = normalizePath(
      file.path(project_root, "stan_cache"),
      mustWork = FALSE
    ))
  }
})

pilot_methods <- c(
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
  # sydni_factor synthesizes a discretized version of the continuous y
  # outcome (engine_factor() requires integer-like/ordinal input); its rows
  # are flagged via outcome_discretized and should be read as directional
  # evidence about the factor engine, not as directly comparable in
  # magnitude to the other methods' continuous-y results. See
  # discretize_to_ordinal() in inferential_safety_framework.R.
  "sydni_factor",
  "sydni_brm"
)

raw_path <- "results/benchmarks/inferential_safety_pilot_raw.csv"
summary_path <- "results/benchmarks/inferential_safety_pilot_summary.csv"

# checkpoint_path makes this resumable: each condition's row is appended to
# raw_path as soon as it finishes, and any condition already present there
# (from a prior, interrupted run) is skipped rather than recomputed. Safe to
# kill this process at any point - rerunning this script picks up wherever
# it left off instead of starting over. Delete raw_path first if you want a
# genuinely fresh run instead of a resume.
pilot_results <- run_inferential_safety_benchmark(
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
  methods = pilot_methods,
  reps = 10,
  p = 12,
  seed = 20260612,
  checkpoint_path = raw_path,
  n_workers = n_workers,
  worker_cores = worker_cores,
  worker_setup = worker_setup
)

pilot_summary <- summarize_inferential_safety(pilot_results)

utils::write.csv(pilot_summary, summary_path, row.names = FALSE)

cat("Wrote raw pilot results to:", raw_path, "\n")
cat("Wrote summarized pilot results to:", summary_path, "\n\n")

headline <- subset(
  pilot_summary,
  method %in% c(
    "SMOTE", "ADASYN", "synthpop_sample", "synthpop_cart",
    "sydni_copula", "sydni_graphical", "sydni_bart", "sydni_factor", "sydni_brm"
  )
)
headline <- headline[order(headline$scenario, headline$aug_ratio, headline$method), ]

if (any(headline$outcome_discretized)) {
  cat("NOTE: rows with method == \"sydni_factor\" (outcome_discretized ==",
      "TRUE) use a discretized version of the continuous y outcome, because",
      "engine_factor() does not accept continuous input. Their rates are",
      "directional evidence about the factor engine, not directly",
      "comparable in magnitude to the continuous-y results for the other",
      "methods below.\n\n")
}

# Two clearly separated tables, matching the benchmark's two distinct
# questions: (1) how bad is the discouraged "pool the data and treat it as a
# bigger-N discovery analysis" usage, tested here as an explicit off-label
# stress test the same way SMOTE/ADASYN are; and (2) how safe/useful is
# SYDNI's actual documented workflow, where the observed analysis stays
# primary and the augmented analysis is only a secondary, cautiously
# interpreted check. These should not be read as one combined verdict.

cat("=== Off-label stress test: pooled data treated as a standalone",
    "discovery analysis (not SYDNI's documented workflow) ===\n\n")
print(
  headline[, c(
    "scenario",
    "n_obs",
    "aug_ratio",
    "method",
    "type1_rate",
    "false_discovery_rate",
    "coverage_95",
    "failed",
    "warning_count",
    "runtime_sec"
  )],
  row.names = FALSE,
  digits = 3
)

cat("\n=== Intended-use safety: observed data primary, augmented data as a",
    "secondary check (SYDNI's documented workflow) ===\n\n")
print(
  headline[, c(
    "scenario",
    "n_obs",
    "aug_ratio",
    "method",
    "false_augmentation_only_sig_rate",
    "power_recovery_rate",
    "mean_effect_ratio_nonzero_vs_obs",
    "mean_std_effect_change_vs_obs"
  )],
  row.names = FALSE,
  digits = 3
)

cat("\n=== Structure retention, shrinkage-aware ===\n\n")
print(
  headline[, c(
    "scenario",
    "n_obs",
    "aug_ratio",
    "method",
    "mean_abs_cor_error",
    "cor_sign_reversal_rate",
    "mean_cor_magnitude_change_nonnull"
  )],
  row.names = FALSE,
  digits = 3
)

