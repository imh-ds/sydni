# Study 3 runner: empirical illustration on the real `grit` dataset. Single
# run per method (no reps, no grid) - see case_study_grit_framework.R for
# the design rationale. Resumable: methods already present in run_summary.csv
# are skipped, so this can be re-invoked (e.g. after deferring brm to a
# quieter time) without redoing completed methods.

user_lib <- Sys.getenv("R_LIBS_USER")
local_lib <- "local_r_lib"
.libPaths(c(user_lib, local_lib, .libPaths()))

invisible(lapply(list.files("R", full.names = TRUE), source))
source("inst/benchmarks/case_study_grit_framework.R")

stan_cache_dir <- normalizePath(file.path(getwd(), "stan_cache"), mustWork = FALSE)
dir.create(stan_cache_dir, recursive = TRUE, showWarnings = FALSE)
options(cmdstanr_write_stan_file_dir = stan_cache_dir)

out_dir <- "results/case_studies/grit_illustration"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_path <- file.path(out_dir, "run_summary.csv")
coef_path <- file.path(out_dir, "downstream_coefficients.csv")
structure_path <- file.path(out_dir, "structure_metrics.csv")

all_methods <- c(
  "observed_only", "bootstrap", "noise_bootstrap", "SMOTE", "ADASYN",
  "synthpop_sample", "synthpop_cart",
  "sydni_copula_flat", "sydni_copula_ma",
  "sydni_graphical_flat", "sydni_graphical_ma",
  "sydni_bart_flat", "sydni_bart_ma",
  "sydni_factor_flat", "sydni_factor_ma",
  "sydni_brm_flat", "sydni_brm_ma"
)

# Restrict to a subset via an environment variable, e.g. to defer the two
# brm methods to a later, quieter run:
#   GRIT_CASE_STUDY_METHODS=observed_only,bootstrap,... Rscript ...
methods_env <- Sys.getenv("GRIT_CASE_STUDY_METHODS", unset = NA)
methods <- if (!is.na(methods_env) && nzchar(methods_env)) {
  strsplit(methods_env, ",")[[1]]
} else {
  all_methods
}

done_methods <- character()
if (file.exists(summary_path)) {
  existing <- utils::read.csv(summary_path, stringsAsFactors = FALSE)
  done_methods <- existing$method[!existing$failed]
}
methods <- setdiff(methods, done_methods)

sim <- prepare_grit_case_study_data()
cat("n_obs (complete cases):", nrow(sim$data), "\n")
cat("class split:", paste(names(table(sim$data$class)), table(sim$data$class), collapse = ", "), "\n")
cat("methods to run:", paste(methods, collapse = ", "), "\n\n")

for (m in methods) {
  cat(sprintf("[%s] running...\n", m))
  res <- run_one_grit_case_study(m, sim, aug_ratio = 1, seed = 20260822, cores = 2)
  cat(sprintf("  done (%.2fs, failed=%s)\n", res$summary$runtime_sec, res$summary$failed))

  append_row <- function(df, path) {
    file_exists <- file.exists(path)
    utils::write.table(df, file = path, sep = ",", row.names = FALSE,
                        col.names = !file_exists, append = file_exists, qmethod = "double")
  }
  append_row(res$summary, summary_path)
  if (!is.null(res$coef_table)) append_row(res$coef_table, coef_path)
  if (!is.null(res$structure)) append_row(res$structure, structure_path)
}

cat("\nAll requested methods processed. Current contents of", summary_path, ":\n")
if (file.exists(summary_path)) {
  print(utils::read.csv(summary_path, stringsAsFactors = FALSE)[, c("method", "runtime_sec", "failed")])
}
