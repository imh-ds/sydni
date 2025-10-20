---
editor_options: 
  markdown: 
    wrap: 72
---

# sydni

`sydni` is an R package for running *SYDNI: Synthetic Data Augmentation
through Noise Inference*. SYDNI is a framework for generating "muted
look-alike" synthetic datasets by fitting a model to the observed data
with controlled noise-supported shrinkage. Its goal is to support
robustness and sensitivity checks, particularly in small or imbalanced
samples, while guarding against inflated Type I error rates (false
positives).

## Use-case

Sample size limitations are common in behavioral and business sciences
for various reasons, including (but not limited to):

-   High costs (e.g., participant compensation, recruitment)

-   Natural imbalances (e.g., gender ratios in human subject pools)

-   Niche or underrepresented populations (e.g., rare diseases)

These constraints limit researchers' ability to generalize findings or
run analyses that require large or balanced datasets (e.g., clustering,
latent models). Such limitations more commonly disproportionately affect
researchers from underfunded institutions or marginalized communities
who lack the resources to overcome the practical obstacles.

**SYDNI provides a principled way to augment data** using noisy
synthetic samples that:

-   Mimic observed structure without overfitting idiosyncrasies

-   Penalize spurious signal

-   Allow for more confident robustness and sensitivity analyses

This makes SYDNI especially useful when working with small samples,
biased groups, or constrained study designs.

## Installation

Install `sydni` from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("imh-ds/sydni")
```

For optimal performance using `cmdstanr` as the backend instead of
`rstan`, you may also need to install `cmdstanr` (for details, see the
CmdStanR website: <https://mc-stan.org/cmdstanr/>):

``` r
remotes::install_github("stan-dev/cmdstanr")

# Ensure CmdStan is also installed
cmdstanr::install_cmdstan()
```

## Quick Start

### Engines

`sydni` currently supports two 'engines' to generate synthetic data:

-   **Bayesian Regression Models (BRM):** Fits a BRM to the observed
    data and applies shrinkage via priors. A noisy resampled dataset is
    used as the 'feeder' data into generating the synthetic data.

-   **Psychometric Copula:** Fits a copula to the noisy data with a
    shrunken correlation matrix.

### Generating Synthetic Data

To generate synthetic data, specify:

-   `data`: The observed data to fit the model to

-   `variables`: The variables to synthesize

-   `n_aug`: The sample size of the desired synthetic dataset

-   `engine`: The model to synthesize data

``` r
library(tidyverse)

# Get the observed dataset
grit <- data(grit) |>
  filter(if_all(everything(), ~ !is.na(.)))

# Run SYDNI
brm_data <- sydni(
  data = grit,
  variables = c(sprintf("o%01d", seq(10)),
                sprintf("e%01d", seq(10))),
  n_aug = 100,
  engine = "brm"
)
```

### Generating Synthetic Data with Covariates

If there are any covariates that should be taken into consideration when
generating synthetic data, these can be specified using the `covariates`
argument within the `sydni()` function, like so:

``` r
# Run SYDNI with covariates
brm_data <- sydni(
  data = grit,
  variables = c(sprintf("o%01d", seq(10)),
                sprintf("e%01d", seq(10))),
  covariates = c("gs1", "gs2"),
  n_aug = 100,
  engine = "copula"
)
```

Covariates do not have any noise added during the synthetic data
generation process and are only resampled from the observed data. These
covariates are appended to the synthetic data and returned with them in
the final output. This is particularly important in cases where
covariates are categorical, such as gender or group membership, where
noise injection can fundamentally change meaningful data
characteristics.

### Adjusting Noise Injection

To mitigate overfitting to idiosyncratic dependence patterns in small
samples, SYDNI injects controlled noise into the resampled observed data
before fitting the model. This aims to help the synthetic data
generalize better and reflect uncertainty. You can control the amount
and type of noise using the arguments below:

-   `method`: One of `"discrete"`, `"truncnormal"`, `"gaussian"`, or
    `"random"` determines the shape of the noise distribution.

-   `c`: Controls the scale of the noise relative to the standard
    deviation of the variable. Used in most noise types. If `"auto"`
    (default), it’s computed based on either p or f below.

-   `p`: For `"discrete"` noise, controls the probability that a value
    changes by ±1. Higher values of `p` result in more perturbation.
    Only used when `c = "auto"`.

-   `f`: For continuous noise (e.g., `"gaussian"`, `"truncnormal"`),
    this sets the target ratio of noise variance to data variance. For
    example, `f = 0.25` injects noise with one-quarter the variance of
    the original variable. Only used when `c = "auto"`.

-   `r`: When using `method = "random"`, this sets the maximum noise
    magnitude added (drawn from `Uniform(-r, r)`), independent of the
    data’s scale.

Inject mild discrete noise:

``` r
sydni(..., engine = "copula", method = "discrete", p = 0.2)
```

Use more aggressive continuous noise:

``` r
sydni(..., engine = "copula", method = "gaussian", f = 0.5) 
```

Disable auto-scaling and specify exact scale:

``` r
sydni(..., engine = "copula", method = "truncnormal", c = 0.1) 
```

Add unscaled uniform noise:

``` r
sydni(..., engine = "copula", method = "random", r = 1)
```

### Adjusting Shrinkage Magnitude (for BRM only)

When running SYDNI with `engine = "brm"`, the effect size estimates are
regularized by placing normal priors centered at 0 on all regression
coefficients. The amount of shrinkage is controlled by two key
arguments: `k` and `base_sd`.

-   `base_sd` controls the starting width of the prior when no synthetic
    data is used.

-   `k` controls how aggressively the prior contracts as more synthetic
    data is introduced.

The prior standard deviation is computed as:

``` r
SD = base_sd / (1 + n_syn / n_obs)^k
```

This means:

-   Higher values of `k` leads to more aggressive shrinkage as synthetic
    data increases.

-   Lower values of `k` (e.g., 0.5–1) leads to more conservative, slower
    shrinkage.

-   Higher `base_sd` allows more flexibility; lower values pull
    estimates tighter around zero.

``` r
# More aggressive shrinkage (e.g., discourage spurious effects)
sydni(..., engine = "brm", k = 3)

# More permissive shrinkage (e.g., allow moderate effects to survive)
sydni(..., engine = "brm", k = 0.5)
```

You can also raise or lower `base_sd` to control the overall scale of
prior uncertainty. For example, setting `base_sd = 5` will produce
tighter priors than the default of 10, even when `k = 1`.

## Citation

If you use `sydni` in published work, please cite:

> Im, H. (2025). sydni: Synthetic Data Augmentation through Noise
> Inference. GitHub repository. <https://github.com/imh-ds/sydni>.
