# SYDNI copula line-by-line real-data example
#
# This script intentionally avoids the benchmark wrapper. It is meant to be run
# interactively, section by section, in a standard R session.

# -------------------------------------------------------------------------
# 0. Libraries and package source
# -------------------------------------------------------------------------

.libPaths(c(
  Sys.getenv("R_LIBS_USER"),
  "local_r_lib",
  .libPaths()
))

library(dplyr)

source_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in source_files) {
  source(file, local = globalenv())
}

set.seed(20260613)


# -------------------------------------------------------------------------
# 1. Helper functions
# -------------------------------------------------------------------------

starts_with_name <- function(nms, prefix) {
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

tidy_lm <- function(fit) {
  tab <- as.data.frame(summary(fit)$coefficients)
  tab$term <- rownames(tab)
  names(tab)[1:4] <- c("estimate", "std_error", "statistic", "p_value")
  rownames(tab) <- NULL
  tab
}


# -------------------------------------------------------------------------
# 2. Load and prepare Experiment 01
# -------------------------------------------------------------------------

exp1_raw <- read.csv(
  "tests/3_Experiments/Experiment_01/1_Data/cc_data_s2_data.csv",
  check.names = FALSE
)

exp1_data <- exp1_raw |>
  dplyr::filter(
    attention == 1,
    Late == 0,
    gen != 3
  )

exp1_keep <- c(
  "ResponseId", "age", "gen", "DAD", "MOM",
  starts_with_name(names(exp1_data), "mss"),
  starts_with_name(names(exp1_data), "dcs"),
  starts_with_name(names(exp1_data), "hssb")
)

exp1_keep <- exp1_keep[!grepl("DO$", exp1_keep)]

exp1_data <- exp1_data[, exp1_keep, drop = FALSE]
names(exp1_data) <- tolower(names(exp1_data))

exp1_data$responseid <- sprintf("id_%03d", seq_len(nrow(exp1_data)))
exp1_data <- median_impute(exp1_data)


# -------------------------------------------------------------------------
# 3. Define variables, subgroup, and analysis model
# -------------------------------------------------------------------------

exp1_variables <- c(
  paste0("mss_", 1:4),
  paste0("dcs_", 1:4),
  paste0("hssb_", 1:12)
)

exp1_covariates <- c("age", "gen", "dad", "mom")
exp1_id <- "responseid"
exp1_group <- "gen"
exp1_group_value <- 2
exp1_n_aug <- 142

exp1_observed <- exp1_data[exp1_data[[exp1_group]] == exp1_group_value, , drop = FALSE]
exp1_target_items <- exp1_observed[, exp1_variables, drop = FALSE]


# -------------------------------------------------------------------------
# 4. Run SYDNI copula
# -------------------------------------------------------------------------

# Stronger discrete noise is used here because these are Likert-style item
# indicators and the purpose is conservative augmentation.
exp1_copula_fit <- sydni(
  data = exp1_target_items,
  variables = exp1_variables,
  covariates = NULL,
  engine = "copula",
  n_aug = exp1_n_aug,
  noise_method = "discrete",
  c = "auto",
  p = 0.80,
  round_int = TRUE,
  seed = 20260613
)

exp1_synthetic_items <- exp1_copula_fit$data


# -------------------------------------------------------------------------
# 5. Add covariates back to the synthetic rows
# -------------------------------------------------------------------------

set.seed(20260613)

exp1_covariate_draws <- exp1_observed[
  sample(seq_len(nrow(exp1_observed)), size = nrow(exp1_synthetic_items), replace = TRUE),
  setdiff(exp1_covariates, exp1_group),
  drop = FALSE
]

exp1_synthetic <- cbind(exp1_synthetic_items, exp1_covariate_draws)
exp1_synthetic[[exp1_group]] <- exp1_group_value
exp1_synthetic[[exp1_id]] <- sprintf("aug_copula_%03d", seq_len(nrow(exp1_synthetic)))

missing_cols <- setdiff(names(exp1_data), names(exp1_synthetic))
for (col in missing_cols) {
  exp1_synthetic[[col]] <- NA
}

exp1_synthetic <- exp1_synthetic[, names(exp1_data), drop = FALSE]
exp1_synthetic <- median_impute(exp1_synthetic)


# -------------------------------------------------------------------------
# 6. Create composites from observed and synthetic indicator-level data
# -------------------------------------------------------------------------

exp1_observed_scored <- add_exp1_composites(exp1_observed)
exp1_synthetic_scored <- add_exp1_composites(exp1_synthetic)

exp1_pooled_scored <- rbind(exp1_observed_scored, exp1_synthetic_scored)


# -------------------------------------------------------------------------
# 7. Run observed-only and observed-plus-synthetic regressions
# -------------------------------------------------------------------------

exp1_formula <- dissatis ~ material + merit + opportunity + chance + age + dad + mom

exp1_observed_lm <- lm(exp1_formula, data = exp1_observed_scored)
exp1_pooled_lm <- lm(exp1_formula, data = exp1_pooled_scored)

summary(exp1_observed_lm)
summary(exp1_pooled_lm)


# -------------------------------------------------------------------------
# 8. Compare predictor-level coefficients
# -------------------------------------------------------------------------

exp1_observed_tab <- tidy_lm(exp1_observed_lm)
exp1_pooled_tab <- tidy_lm(exp1_pooled_lm)

exp1_comparison <- merge(
  exp1_observed_tab,
  exp1_pooled_tab,
  by = "term",
  suffixes = c("_observed", "_pooled")
)

exp1_comparison <- exp1_comparison[exp1_comparison$term != "(Intercept)", ]

exp1_comparison$beta_change <- (
  exp1_comparison$estimate_pooled -
    exp1_comparison$estimate_observed
)

exp1_comparison$beta_abs_difference <- abs(
  exp1_comparison$estimate_pooled -
    exp1_comparison$estimate_observed
)

exp1_comparison$effect_magnitude_change <- (
  abs(exp1_comparison$estimate_pooled) -
    abs(exp1_comparison$estimate_observed)
)

exp1_comparison$effect_magnitude_inflation <- pmax(
  exp1_comparison$effect_magnitude_change,
  0
)

exp1_comparison$effect_magnitude_shrinkage <- pmax(
  -exp1_comparison$effect_magnitude_change,
  0
)

exp1_comparison$se_ratio <- (
  exp1_comparison$std_error_pooled /
    exp1_comparison$std_error_observed
)

exp1_comparison$augmentation_only_sig <- (
  exp1_comparison$p_value_observed >= 0.05 &
    exp1_comparison$p_value_pooled < 0.05
)

exp1_comparison$sign_flip <- (
  sign(exp1_comparison$estimate_observed) !=
    sign(exp1_comparison$estimate_pooled)
)

exp1_comparison <- exp1_comparison[
  order(exp1_comparison$effect_magnitude_inflation, decreasing = TRUE),
]

print(exp1_comparison)


# -------------------------------------------------------------------------
# 9. Optional: save output
# -------------------------------------------------------------------------

dir.create("results/benchmarks", recursive = TRUE, showWarnings = FALSE)

write.csv(
  exp1_synthetic,
  "results/benchmarks/sydni_copula_exp1_synthetic_rows.csv",
  row.names = FALSE
)

write.csv(
  exp1_comparison,
  "results/benchmarks/sydni_copula_exp1_coefficient_comparison.csv",
  row.names = FALSE
)
