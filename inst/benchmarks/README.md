# SYDNI inferential-safety benchmarks

This folder contains benchmark scaffolds for evaluating SYDNI's intended use:
conservative synthetic augmentation for robustness and sensitivity analysis.

The main benchmark target is **inferential safety under off-label data
augmentation**, not classifier balancing, synthetic training-data utility, or
privacy-preserving microdata release. In particular, comparisons with SMOTE,
ADASYN, GAN-style tabular generators, and `synthpop` should be framed as
off-label augmentation stress tests:

- SMOTE and ADASYN are included because analysts may repurpose them for sample
  augmentation, even though their native goal is classifier balancing.
- `synthpop` is included because analysts may repurpose synthetic microdata
  generators for augmentation, even though its native goal is disclosure-aware
  synthetic data generation.
- GAN-style methods may be added as optional generative-ML comparators when a
  reproducible local setup is available.

## Files

- `simulation_framework.R`: minimal legacy scaffold for quick experiments.
- `inferential_safety_framework.R`: replicated simulation framework for the
  paper benchmark.

## Recommended benchmark object

The primary analysis target should be the **pooled augmented dataset**:

```r
observed_plus_synthetic <- rbind(observed, synthetic)
```

Synthetic-only results should be retained as diagnostics, but the main question
is whether augmentation makes ordinary statistical analyses appear stronger,
more precise, or more significant than the observed data justify.

## Primary metrics

The paper benchmark should prioritize:

- Type I error for true-null coefficients
- false discovery rate across many null predictors
- augmentation-only significance rate
- effect-size inflation relative to truth and observed-only estimates
- standard-error compression
- 95% interval coverage
- sign reversal for true nonzero effects
- moderation false-positive rate for true-null interactions

Secondary metrics should cover correlation preservation, Likert category
fidelity, psychometric structure, subgroup proportions, duplicate rates,
degenerate outputs, warning counts, and runtime.

## Example usage

```r
source("inst/benchmarks/inferential_safety_framework.R")

pilot <- run_inferential_safety_benchmark(
  reps = 50,
  n_obs_values = c(75, 150),
  aug_ratios = c(1),
  methods = c(
    "observed_only",
    "bootstrap",
    "noise_bootstrap",
    "SMOTE",
    "ADASYN",
    "synthpop_sample",
    "sydni_copula",
    "sydni_graphical"
  )
)

summary <- summarize_inferential_safety(pilot)
summary
```

For manuscript-scale benchmarks, use at least 500 replications per condition
after the pilot is stable.
