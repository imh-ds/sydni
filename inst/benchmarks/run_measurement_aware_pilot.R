# Run Simulation Study 2: the measurement-aware (two-stage) calibration
# benchmark. See inst/benchmarks/measurement_aware_framework.R for the
# rationale and design; results/benchmarks/simulation_study_2/ for output.
#
# This is a separate simulation from Simulation Study 1
# (inst/benchmarks/run_inferential_safety_pilot.R / results in
# results/benchmarks/*.csv and FINDINGS.md), which never exercised
# measurement_blocks. Same n_obs/aug_ratio grid and replicate count as
# Study 1's pilot, for comparability.

user_lib <- Sys.getenv("R_LIBS_USER")
local_lib <- "local_r_lib"
.libPaths(c(user_lib, local_lib, .libPaths()))

invisible(lapply(list.files("R", full.names = TRUE), source))
source("inst/benchmarks/measurement_aware_framework.R")

dir.create("results/benchmarks/simulation_study_2", recursive = TRUE, showWarnings = FALSE)

stan_cache_dir <- normalizePath(file.path(getwd(), "stan_cache"), mustWork = FALSE)
dir.create(stan_cache_dir, recursive = TRUE, showWarnings = FALSE)
options(cmdstanr_write_stan_file_dir = stan_cache_dir)

# Parallel dispatch, same rationale as Study 1 (see
# run_inferential_safety_pilot.R): sydni_brm dominates runtime and caps
# itself to 2 cores per call (chains = 2), so concurrent independent
# conditions are close to embarrassingly parallel. n_workers = 4 leaves 4
# of this machine's 12 physical cores free.
n_workers <- 4
worker_cores <- 2
project_root <- getwd()
worker_setup <- local({
  root <- project_root
  function() {
    setwd(root)
    .libPaths(c(
      Sys.getenv("R_LIBS_USER"),
      "local_r_lib",
      .libPaths()
    ))
    invisible(lapply(list.files("R", full.names = TRUE), source))
    source("inst/benchmarks/measurement_aware_framework.R")
    options(cmdstanr_write_stan_file_dir = normalizePath(
      file.path(root, "stan_cache"),
      mustWork = FALSE
    ))
  }
})

# 17 arms: the 7 off-label/incumbent comparators from Study 1, plus each
# of the 5 measurement-aware-capable SYDNI engines run twice - "_flat"
# (no measurement_blocks, identical call shape to every engine call in
# Study 1) and "_ma" (measurement-aware, with measurement_blocks +
# between_block_shrink). The flat/ma pair for a given engine, run on
# identical simulated data with identical seeds, isolates what
# measurement_blocks itself contributes - see
# measurement_aware_framework.R's augment_method_measurement_aware().
pilot_methods <- c(
  "observed_only",
  "bootstrap",
  "noise_bootstrap",
  "SMOTE",
  "ADASYN",
  "synthpop_sample",
  "synthpop_cart",
  "sydni_copula_flat", "sydni_copula_ma",
  "sydni_graphical_flat", "sydni_graphical_ma",
  "sydni_bart_flat", "sydni_bart_ma",
  "sydni_factor_flat", "sydni_factor_ma",
  "sydni_brm_flat", "sydni_brm_ma"
)

raw_path <- "results/benchmarks/simulation_study_2/measurement_aware_pilot_raw.csv"
summary_path <- "results/benchmarks/simulation_study_2/measurement_aware_pilot_summary.csv"

pilot_results <- run_measurement_aware_benchmark(
  n_obs_values = c(75, 150, 300),
  aug_ratios = c(0.5, 1, 2),
  methods = pilot_methods,
  reps = 10,
  seed = 20260821,
  checkpoint_path = raw_path,
  n_workers = n_workers,
  worker_cores = worker_cores,
  worker_setup = worker_setup
)

cat("Wrote raw pilot results to:", raw_path, "\n")
cat("Total rows:", nrow(pilot_results), "\n")
cat("Failed rows:", sum(pilot_results$failed, na.rm = TRUE), "\n\n")

pilot_summary <- aggregate(
  cbind(
    type1_rate, false_discovery_rate, false_augmentation_only_sig_rate,
    power_recovery_rate, within_scale_alpha_abs_change, within_scale_cor_error,
    between_construct_cor_error, between_construct_sign_reversal_rate,
    between_construct_magnitude_change, runtime_sec
  ) ~ method + n_obs + aug_ratio,
  data = pilot_results[pilot_results$failed != TRUE, ],
  FUN = function(x) mean(x, na.rm = TRUE)
)
utils::write.csv(pilot_summary, summary_path, row.names = FALSE)
cat("Wrote summarized pilot results to:", summary_path, "\n")
