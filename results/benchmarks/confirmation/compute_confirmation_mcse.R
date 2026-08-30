# Design-stratified Monte Carlo standard errors (MCSEs) and 95% Monte Carlo
# interval (MCI) half-widths for the confirmation-run headline estimates
# reported in Sections 8.6-8.7 of paper/sydni-jss-methodology-draft.md.
#
# WHY "design-stratified": every headline number in this study is an
# unweighted average of per-replicate metric values across a factorial grid
# of design cells (scenario x n_obs x aug_ratio, or n_obs x aug_ratio for
# Study 2). Pooling the raw SD across every row in that average and dividing
# by sqrt(total rows) conflates two different sources of variation: real,
# by-design heterogeneity across cells (different scenarios/sample sizes are
# *supposed* to have different true rates) and genuine Monte Carlo sampling
# error (repeated replicates of the same cell varying by chance). Only the
# second is relevant to "how precisely do we know this mean." This script
# instead computes, for each design cell, the replicate-level variance and
# divides by that cell's replicate count, then combines cells with equal
# weight (matching how the paper's grand means are computed: an unweighted
# mean over cells with equal replicate counts per cell equals an unweighted
# mean over rows). See Morris, White & Crowther (2019, Stat Med) for the
# general Monte Carlo standard error framework this follows.
#
# Usage: Rscript results/benchmarks/confirmation/compute_confirmation_mcse.R
# Output: results/benchmarks/confirmation/study1_confirmation_mcse.csv
#         results/benchmarks/confirmation/study2_confirmation_mcse.csv
# Both are also printed to the console in the grouping used by the paper's
# tables, so the numbers below can be checked directly against the tables in
# paper/sydni-jss-methodology-draft.md Sections 8.6-8.7.

z975 <- qnorm(0.975)

#' Design-stratified MCSE of an unweighted grand mean of `metric_col` across
#' cells defined by `strata_cols`, restricted by an optional row filter
#' already applied to `data` before calling this function.
#'
#' Assumes (and warns if violated) that every stratum has the same replicate
#' count, which is true of every grid used in this confirmation study - the
#' grand mean is then exactly the unweighted mean of the per-stratum means,
#' and its variance is the average per-stratum replicate variance divided by
#' (n_strata * reps_per_stratum).
stratified_mcse <- function(data, strata_cols, metric_col) {
  data <- data[!is.na(data[[metric_col]]), ]
  cell_id <- interaction(data[strata_cols], drop = TRUE)
  by_cell <- split(data[[metric_col]], cell_id)
  reps <- vapply(by_cell, length, integer(1))
  if (length(unique(reps)) > 1) {
    warning(sprintf(
      "Unbalanced replicate counts across strata (%s): %s",
      paste(strata_cols, collapse = "/"), paste(unique(reps), collapse = ",")
    ))
  }
  cell_means <- vapply(by_cell, mean, numeric(1))
  cell_vars <- vapply(by_cell, function(x) if (length(x) > 1) stats::var(x) else NA_real_, numeric(1))
  h <- length(by_cell)
  grand_mean <- mean(cell_means)
  # Var(grand mean) = (1/H^2) * sum_h(cell_var_h / reps_h); NA cell variances
  # (a stratum with a single successful replicate) are dropped from the sum
  # with a warning rather than silently propagating NA.
  valid <- !is.na(cell_vars)
  if (any(!valid)) {
    warning(sprintf(
      "%d of %d strata have <2 successful replicates for %s; excluded from MCSE sum (mean still uses all strata).",
      sum(!valid), h, metric_col
    ))
  }
  var_grand <- sum(cell_vars[valid] / reps[valid]) / h^2
  mcse <- sqrt(var_grand)
  data.frame(
    metric = metric_col,
    mean = grand_mean,
    mcse = mcse,
    ci95_halfwidth = z975 * mcse,
    n_strata = h,
    reps_per_stratum = paste(sort(unique(reps)), collapse = ","),
    n_rows = sum(reps)
  )
}

#' Convenience wrapper: run stratified_mcse() for one method at a time and
#' rbind the results, adding a `method` column.
by_method_mcse <- function(data, method_col, strata_cols, metric_cols) {
  methods <- sort(unique(data[[method_col]]))
  out <- do.call(rbind, lapply(methods, function(m) {
    dm <- data[data[[method_col]] == m, ]
    rows <- do.call(rbind, lapply(metric_cols, function(mc) stratified_mcse(dm, strata_cols, mc)))
    cbind(method = m, rows)
  }))
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Study 1 (Section 8.6)
# ---------------------------------------------------------------------------

s1_nonbrm <- read.csv("results/benchmarks/confirmation/study1_nonbrm_confirmation_raw.csv", stringsAsFactors = FALSE)
s1_brm <- read.csv("results/benchmarks/confirmation/study1_brm_confirmation_raw.csv", stringsAsFactors = FALSE)
s1_nonbrm <- s1_nonbrm[s1_nonbrm$failed != TRUE, ]
s1_brm <- s1_brm[s1_brm$failed != TRUE, ]

matched_scenarios <- c("null_only", "sparse_effects", "correlated_predictors", "moderation_imbalance", "small_cell")
matched_n_obs <- c(75, 150)
matched_aug <- c(1, 2)
s1_strata <- c("scenario", "n_obs", "aug_ratio")

s1_nonbrm_matched <- s1_nonbrm[s1_nonbrm$scenario %in% matched_scenarios &
                                  s1_nonbrm$n_obs %in% matched_n_obs &
                                  s1_nonbrm$aug_ratio %in% matched_aug, ]

# Engine-tradeoff table (Section 8.6): 4 continuous-outcome engines + sydni_brm,
# on the matched grid, for the 3 metrics printed in that table.
tradeoff_engines <- c("sydni_graphical", "sydni_copula", "sydni_bart")
tradeoff_metrics <- c("type1_rate", "power_recovery_rate", "null_edge_spurious_rate")

tradeoff_nonbrm <- by_method_mcse(
  s1_nonbrm_matched[s1_nonbrm_matched$method %in% tradeoff_engines, ],
  "method", s1_strata, tradeoff_metrics
)
tradeoff_brm <- by_method_mcse(
  s1_brm[s1_brm$method == "sydni_brm", ],
  "method", s1_strata, tradeoff_metrics
)
engine_tradeoff_table <- rbind(tradeoff_nonbrm, tradeoff_brm)

cat("\n=== Section 8.6 engine-tradeoff table (matched grid) ===\n")
print(engine_tradeoff_table[order(engine_tradeoff_table$metric, engine_tradeoff_table$method),
                             c("method", "metric", "mean", "mcse", "ci95_halfwidth", "n_rows")],
      row.names = FALSE, digits = 4)

# sydni_factor, reported separately: same 3 metrics on the matched grid, plus
# false_augmentation_only_sig_rate on the FULL (non-matched) grid.
factor_matched <- by_method_mcse(
  s1_nonbrm_matched[s1_nonbrm_matched$method == "sydni_factor", ],
  "method", s1_strata, tradeoff_metrics
)
factor_full_strata <- c("scenario", "n_obs", "aug_ratio")
factor_full <- by_method_mcse(
  s1_nonbrm[s1_nonbrm$method == "sydni_factor", ],
  "method", factor_full_strata, "false_augmentation_only_sig_rate"
)
factor_table <- rbind(factor_matched, factor_full)

cat("\n=== sydni_factor's own numbers (Section 8.6) ===\n")
print(factor_table[, c("method", "metric", "mean", "mcse", "ci95_halfwidth", "n_rows")],
      row.names = FALSE, digits = 4)

# Intended-use safety paragraph: false_augmentation_only_sig_rate, full grid,
# 4 continuous-outcome engines.
safety_engines <- c("sydni_graphical", "sydni_copula", "sydni_bart")
safety_table <- by_method_mcse(
  s1_nonbrm[s1_nonbrm$method %in% safety_engines, ],
  "method", factor_full_strata, "false_augmentation_only_sig_rate"
)

cat("\n=== Intended-use safety: false_augmentation_only_sig_rate (full grid) ===\n")
print(safety_table[, c("method", "mean", "mcse", "ci95_halfwidth", "n_rows")],
      row.names = FALSE, digits = 4)

# Off-label stress test: pooled Type I error, SYDNI+synthpop_cart side vs
# off-label side, full grid. Strata here are scenario x n_obs x aug_ratio x
# method, since methods on the same side of the comparison have genuinely
# different true rates by design and shouldn't be pooled into "noise."
sydni_side <- c("sydni_copula", "sydni_graphical", "sydni_bart", "sydni_factor", "synthpop_cart")
offlabel_side <- c("bootstrap", "noise_bootstrap", "SMOTE", "ADASYN")
offlabel_strata <- c("scenario", "n_obs", "aug_ratio", "method")

sydni_side_mcse <- stratified_mcse(s1_nonbrm[s1_nonbrm$method %in% sydni_side, ], offlabel_strata, "type1_rate")
offlabel_side_mcse <- stratified_mcse(s1_nonbrm[s1_nonbrm$method %in% offlabel_side, ], offlabel_strata, "type1_rate")
offlabel_comparison <- rbind(
  cbind(side = "sydni_synthpop_cart", sydni_side_mcse),
  cbind(side = "offlabel_comparators", offlabel_side_mcse)
)

cat("\n=== Off-label stress test: pooled Type I error, method x scenario x n_obs x aug_ratio strata ===\n")
print(offlabel_comparison[, c("side", "mean", "mcse", "ci95_halfwidth", "n_rows")],
      row.names = FALSE, digits = 4)

# sydni_brm's own "every primary metric under 1pp" claim (Section 8.6 intro):
# the four metrics that appear in its rows across the paper (the three
# engine-tradeoff metrics above, computed already, plus
# false_augmentation_only_sig_rate on its own narrower grid).
brm_fasr <- stratified_mcse(s1_brm[s1_brm$method == "sydni_brm", ], s1_strata, "false_augmentation_only_sig_rate")
cat("\n=== sydni_brm's false_augmentation_only_sig_rate (its own narrow grid) ===\n")
print(brm_fasr[, c("mean", "mcse", "ci95_halfwidth", "n_rows")], row.names = FALSE, digits = 4)

s1_out <- rbind(
  cbind(table = "engine_tradeoff", engine_tradeoff_table[, c("method", "metric", "mean", "mcse", "ci95_halfwidth", "n_strata", "reps_per_stratum", "n_rows")]),
  cbind(table = "sydni_factor_own", factor_table[, c("method", "metric", "mean", "mcse", "ci95_halfwidth", "n_strata", "reps_per_stratum", "n_rows")]),
  cbind(table = "intended_use_safety", safety_table[, c("method", "metric", "mean", "mcse", "ci95_halfwidth", "n_strata", "reps_per_stratum", "n_rows")]),
  cbind(table = "offlabel_stress_test", data.frame(method = offlabel_comparison$side, offlabel_comparison[, c("metric", "mean", "mcse", "ci95_halfwidth", "n_strata", "reps_per_stratum", "n_rows")])),
  cbind(table = "brm_fasr_own_grid", data.frame(method = "sydni_brm", brm_fasr[, c("metric", "mean", "mcse", "ci95_halfwidth", "n_strata", "reps_per_stratum", "n_rows")]))
)
write.csv(s1_out, "results/benchmarks/confirmation/study1_confirmation_mcse.csv", row.names = FALSE)
cat("\nWrote results/benchmarks/confirmation/study1_confirmation_mcse.csv (", nrow(s1_out), "rows)\n")

# ---------------------------------------------------------------------------
# Study 2 (Section 8.7, Scenario 3 ma-engine table)
# ---------------------------------------------------------------------------

s2_nonbrm <- read.csv("results/benchmarks/confirmation/study2_nonbrm_confirmation_raw.csv", stringsAsFactors = FALSE)
s2_brm <- read.csv("results/benchmarks/confirmation/study2_brm_confirmation_raw.csv", stringsAsFactors = FALSE)
s2_nonbrm <- s2_nonbrm[s2_nonbrm$failed != TRUE, ]
s2_brm <- s2_brm[s2_brm$failed != TRUE, ]

s2_matched_n_obs <- c(75, 150)
s2_matched_aug <- c(1, 2)
s2_strata <- c("n_obs", "aug_ratio")

ma_engines <- c("sydni_factor_ma", "sydni_graphical_ma", "sydni_copula_ma", "sydni_bart_ma")
ma_metrics <- c("type1_rate", "false_discovery_rate", "power_recovery_rate",
                 "within_scale_cor_error", "between_construct_sign_reversal_rate", "runtime_sec")
ma_metrics <- ma_metrics[ma_metrics %in% names(s2_nonbrm)]

s2_nonbrm_matched <- s2_nonbrm[s2_nonbrm$n_obs %in% s2_matched_n_obs & s2_nonbrm$aug_ratio %in% s2_matched_aug, ]

scenario3_nonbrm <- by_method_mcse(
  s2_nonbrm_matched[s2_nonbrm_matched$method %in% ma_engines, ],
  "method", s2_strata, ma_metrics
)
scenario3_brm <- by_method_mcse(
  s2_brm[s2_brm$method == "sydni_brm_ma", ],
  "method", s2_strata, ma_metrics
)
scenario3_table <- rbind(scenario3_nonbrm, scenario3_brm)

cat("\n=== Section 8.7 Scenario 3 ma-engine table (matched grid) ===\n")
print(scenario3_table[order(scenario3_table$metric, scenario3_table$method),
                       c("method", "metric", "mean", "mcse", "ci95_halfwidth", "n_rows")],
      row.names = FALSE, digits = 4)

s2_out <- cbind(table = "scenario3_ma_engines", scenario3_table[, c("method", "metric", "mean", "mcse", "ci95_halfwidth", "n_strata", "reps_per_stratum", "n_rows")])
write.csv(s2_out, "results/benchmarks/confirmation/study2_confirmation_mcse.csv", row.names = FALSE)
cat("\nWrote results/benchmarks/confirmation/study2_confirmation_mcse.csv (", nrow(s2_out), "rows)\n")
