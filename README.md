---
editor_options:
  markdown:
    wrap: 72
---

# sydni

`sydni` is an R package for **synthetic data augmentation with
noise-informed shrinkage**. It is designed for **robustness and
sensitivity checks** in settings where analyses may be affected by small
samples, subgroup imbalance, or modest perturbations of observed
structure.

## Scope and intended use

The package is intended to support:

- A tool for generating perturbed synthetic datasets from an observed
  sample
- Sensitivity analyses for evaluating whether downstream conclusions are
  affected by sampling noise and distributional assumptions
- Five synthesis engines: conditional Bayesian regression, Gaussian
  copula with empirical margins, nonlinear BART, latent-factor synthesis,
  and sparse graphical synthesis

The package should not be interpreted as:

- A standalone guarantee of privacy or disclosure protection
- A general solution for recovering population truth from undersized
  samples
- A substitute for domain-specific validation of synthetic data quality

## Installation

Install `sydni` from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("imh-ds/sydni")
```

The BRM engine is optional. Install `brms` only if you want to use
`engine = "brm"`:

```r
install.packages("brms")
```

If you want to use the Bayesian regression engine with the `cmdstanr`
backend:

```r
remotes::install_github("stan-dev/cmdstanr")
cmdstanr::install_cmdstan()
```

Depending on your local Stan and `brms` setup, the BRM engine may also
need compiled dependencies that are installed transitively through
`brms` / `rstan`. If BRM setup does not complete automatically on your
machine, common missing packages include:

```r
install.packages(c("BH", "RcppEigen"))
```

These are not direct `sydni` dependencies, but they are sometimes needed
to finish configuring the BRM toolchain on Windows or on fresh R
installations.

Additional engine package requirements:

```r
install.packages(c("dbarts", "psych", "huge"))
```

- `dbarts` is required for `engine = "bart"`
- `psych` is required for `engine = "factor"`
- `huge` is required for `engine = "graphical"`

## Engines

### Bayesian regression (`engine = "brm"`)

- Builds a noisy bootstrap feeder dataset
- Fits one conditional Bayesian model per synthesized variable on that
  feeder dataset
- Uses shrinkage priors that tighten as the synthetic sample grows
- Best when conditional structure matters and you want model-based
  synthesis
- Provides the most explicit regression-based representation of
  conditional relationships, with support for mixed outcomes and
  interpretable model structure

### Gaussian copula (`engine = "copula"`)

- Builds a noisy feeder dataset within covariate-defined groups
- Estimates subgroup copulas with shrinkage toward a global correlation
  target
- Reconstructs each margin empirically instead of assuming normal output
- Best when dependence preservation is more important than a fully
  specified regression model for each variable
- Emphasizes broad joint dependence structure without fitting a separate
  conditional model for each synthesized variable

### Bayesian additive regression trees (`engine = "bart"`)

- Builds a noisy bootstrap feeder dataset
- Fits one regularized BART model per synthesized variable
- Handles covariates natively like the BRM engine
- Treats ordinal-like synthesized targets through a numeric approximation
  path and reconstructs them back to the observed support with a warning
- Best when you need nonlinear conditional structure
- Extends the conditional-modeling approach to nonlinear relationships
  and interactions

### Latent factor (`engine = "factor"`)

- Builds noisy feeder data within covariate-defined groups
- Fits subgroup factor models for ordinal or integer-like item batteries
- Shrinks factor loadings as augmentation grows
- Best for psychometric indicators and latent-structure robustness checks
- Models latent psychometric structure directly, rather than relying only
  on observed-variable associations

### Sparse graphical (`engine = "graphical"`)

- Builds noisy feeder data within covariate-defined groups
- Fits sparse Gaussian graphical models in a latent Gaussian space
- Reconstructs each margin empirically
- Best when direct conditional-dependence structure is the scientific
  target
- Targets sparse direct-dependence structure rather than dense association
  or latent-factor structure

## Choosing an engine

Use `engine = "brm"` when:

- you want a general-purpose model-based engine
- outcomes are mixed or include ordinal variables that should be handled
  more explicitly than BART's approximation path
- interpretability of conditional structure matters

Use `engine = "copula"` when:

- your main goal is preserving overall dependence structure
- subgroup-aware joint synthesis matters more than explicit regression
  formulas
- you want a relatively assumption-light engine compared with BRM

Use `engine = "bart"` when:

- you suspect important nonlinearities or interactions
- a flexible conditional model is preferred while retaining native
  covariate conditioning
- your outcomes are mostly continuous or binary, or ordinal variables
  are acceptable under an approximation path
- you want tree-based conditional synthesis without switching to a
  dependence-first engine

Use `engine = "factor"` when:

- the synthesized variables are psychometric items or item batteries
- the scientific target is latent structure rather than only observed
  correlations
- you want a robustness check tied to common-factor style structure
- you are working in a measurement setting where latent constructs are a
  more natural model than direct observed-variable links

Use `engine = "graphical"` when:

- the scientific target is direct conditional dependence
- you care about network-like or sparse structure among variables
- a sparse dependency model is preferred to a dense copula correlation
  model
- you want to preserve a sparse dependency skeleton rather than a full
  dense association pattern

## Engine tradeoffs

- `brm` and `bart` are conditional engines. Both model each target from
  the others and handle covariates natively; `brm` emphasizes explicit
  regression structure, whereas `bart` emphasizes nonlinear relationships
  and interactions.
- `copula` and `graphical` are dependence-oriented engines. Both focus on
  joint structure and use subgrouping when covariates are supplied;
  `copula` targets broad association patterns, whereas `graphical`
  targets sparse direct-dependence structure.
- `factor` is the most domain-specific engine. It is intended for item
  batteries and latent constructs, where common-factor structure is more
  central than either regression-based conditioning or dense observed
  correlations.

## Quick start

```r
library(dplyr)
library(sydni)

data(grit, package = "sydni")

grit_complete <- grit |>
  filter(if_all(everything(), ~ !is.na(.)))

result <- sydni(
  data = grit_complete,
  variables = c(sprintf("o%01d", seq(10)), sprintf("e%01d", seq(10))),
  covariates = c("gs1", "gs2"),
  n_aug = 100,
  engine = "copula",
  seed = 123
)

head(result$data)
```

The returned object has three components:

- `models`: fitted model objects or engine-specific model summaries
- `data`: the synthetic dataset
- `meta`: lightweight metadata with generation settings and an input
  summary, not the full original dataset

## Automatic shrinkage and group balancing

When measurement blocks are supplied, omitted shrinkage is calculated from
the effective augmentation:

```r
lambda <- 1 / (1 + n_aug / nrow(data))
```

Here, `n_aug` always means the number of synthetic rows returned for analysis.
If SYDNI generates a larger candidate pool to retain a target group, that
candidate count does not make shrinkage more aggressive. Use
`between_block_shrink = "auto"` to request the same rule explicitly when no
measurement blocks are supplied; synthesized variables are then treated as
singleton structural families.

For an observed sample containing 80 women and 20 men, this call returns the
60 synthetic men required to match the largest observed group:

```r
balanced <- sydni(
  data = dat,
  variables = item_vars,
  covariates = "gender",
  engine = "copula",
  balance_by = "gender",
  balance_value = "man",
  between_block_shrink = "auto",
  seed = 123
)
```

For a manually chosen final augmentation size, use `n_aug`, `retain_by`, and
`retain_value`. Metadata records the candidate count, retained count,
effective augmentation ratio, shrinkage source, and post-retention
calibration diagnostics.

## Noise controls

SYDNI perturbs the feeder data before fitting or sampling.

- `noise_method = "discrete"`:
  for integer-like data, adds `-1`, `0`, or `+1`; for continuous data,
  adds uniform noise scaled by the variable SD
- `noise_method = "gaussian"`:
  adds Gaussian noise scaled by the variable SD
- `noise_method = "truncnormal"`:
  adds truncated Gaussian noise scaled by the variable SD
- `noise_method = "random"`:
  adds uniform noise in `[-r, r]`

Key tuning parameters:

- `p`: target probability of a non-zero discrete perturbation when
  `c = "auto"`
- `f`: target variance ratio for continuous noise when `c = "auto"`
- `r`: direct scale for `noise_method = "random"`
- `noise_bounds`: optional lower/upper clamps for feeder values

Examples:

```r
sydni(..., engine = "copula", noise_method = "discrete", p = 0.2)
sydni(..., engine = "copula", noise_method = "gaussian", f = 0.25)
sydni(..., engine = "copula", noise_method = "truncnormal", c = 0.1)
sydni(..., engine = "copula", noise_method = "random", r = 1)
```

## BRM shrinkage controls

For `engine = "brm"`, slope priors shrink according to:

```r
base_sd / (1 + n_syn / n_obs)^k
```

- Higher `k` means more aggressive shrinkage as synthetic sample size
  increases
- Lower `base_sd` means tighter priors overall
- `sd_scale` controls the intercept prior width relative to the observed
  outcome SD

Example:

```r
sydni(..., engine = "brm", k = 3, base_sd = 5)
```

BRM-specific practical notes:

- `backend = "cmdstanr"` is generally more stable than `backend =
  "rstan"`, especially on Windows.
- The first BRM run can take much longer than later runs because Stan
  model compilation is cached.
- Fast approximate algorithms such as `"meanfield"` may emit warnings
  about unreliable resampling. These warnings usually reflect limitations
  of the approximation rather than an implementation error.

## High-risk use cases

Use extra caution when:

- subgroup cells are extremely small
- variables are heavily multimodal or structurally bounded
- ordinal Likert data are interpreted as fully metric without checking
  the output
- covariate imbalance is severe
- privacy or disclosure protection is a primary goal

## Validation guidance

SYDNI output should be validated before interpretation. At minimum,
compare:

- marginal distributions
- subgroup proportions
- correlations or dependence structure
- downstream model conclusions

## Interpretation and reporting guidelines

SYDNI should be interpreted as a **secondary robustness framework**, not
as additional empirical evidence. The augmented data are generated from
the observed sample and are meant to test whether conclusions appear
stable under cautious, structure-preserving perturbation.

Recommended interpretation rules:

- Treat the original-sample analysis as primary.
- Treat SYDNI analyses as secondary sensitivity or robustness analyses.
- If an effect is present in both the original and augmented analyses,
  that is more reassuring than either analysis alone.
- If an effect is absent in the original data but appears only after
  augmentation, interpret the augmented result as fragile or suggestive,
  not confirmatory.
- Do not present augmentation-only significance as evidence of a latent
  effect that the original data failed to show.

Recommended reporting practices:

- Report the raw-data analysis first.
- Report the SYDNI analysis alongside it, not in place of it.
- Report the full SYDNI settings used:
  `engine`, `n_aug`, `noise_method`, `p`/`f`/`c`, `predictor_method`,
  `predictor_k`, `base_sd`, `k`, `algorithm`, `backend`, and `seed`.
- Report whether conclusions are stable across reasonable alternative
  SYDNI settings rather than emphasizing a single favorable run.
- Describe SYDNI as a robustness or sensitivity procedure, not as a way
  to increase evidential sample size.

Suggested paper language:

> SYDNI results were used as a secondary robustness analysis alongside
> the original-sample results. Synthetic augmentation was not treated as
> additional empirical evidence, and effects that appeared only after
> augmentation were interpreted cautiously.

The package now includes a starter benchmark framework under
`inst/benchmarks/` and a methodology vignette under `vignettes/`.

## Citation

If you use `sydni` in published work, please cite:

> Im, H. (2025). *sydni: Synthetic data augmentation with
> noise-informed shrinkage*. GitHub repository.
> <https://github.com/imh-ds/sydni>
