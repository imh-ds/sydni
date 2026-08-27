# SYDNI: Synthetic Data Augmentation with Noise-Informed Shrinkage for Robustness Analysis in R

Hohjin Im  
University of California, Irvine

## Abstract

Synthetic data methods are often motivated by missing-data handling,
privacy-preserving release, simulation, or generative modeling. In applied
statistical workflows, however, analysts frequently face a narrower but
important question: whether an empirical conclusion is stable when the observed
sample is cautiously perturbed and extended in a way that approximately
preserves relevant marginal, dependence, subgroup, or latent structure.
`sydni` implements Synthetic Data Augmentation with Noise-Informed Shrinkage
(SYDNI), an R package for conservative synthetic augmentation intended for
robustness and sensitivity analysis. SYDNI constructs noisy bootstrap feeder
data, fits one of several swappable synthesis engines, and applies
augmentation-linked shrinkage so that synthetic rows are not interpreted as
ordinary additional empirical evidence. The package currently provides five
engines: conditional Bayesian regression, Gaussian copula synthesis with
empirical margins, Bayesian additive regression trees, latent factor synthesis,
and sparse graphical synthesis. All five engines optionally support
measurement-aware calibration, which protects within-scale item structure
while shrinking between-construct relations, and retention-aware
augmentation, which generates a larger internal candidate pool to return a
targeted, group-balanced synthetic subset. This draft describes the
methodological premise, the general workflow, engine-specific algorithms,
measurement-aware calibration, retention-aware augmentation, validation
recommendations, and a benchmark design suitable for assessing marginal
fidelity, dependence preservation, downstream stability, false-positive
behavior, and computational cost, together with pilot-scale benchmark
results from two simulation studies - a manifest-variable inferential-safety
benchmark, and a companion benchmark isolating the measurement-aware
calibration mechanism on a latent-variable item-battery design - and a
single-run empirical illustration on a real psychometric dataset that
corroborates both studies' central findings.

**Keywords:** synthetic data, robustness analysis, sensitivity analysis,
Bayesian regression, copula, BART, factor analysis, graphical models,
measurement invariance, R.

## 1. Introduction

Many applied analyses are conducted in settings where the observed data are
small, subgroup cells are imbalanced, or complex downstream models are
sensitive to modest perturbations of the empirical distribution. In such
settings, the central inferential concern is often not whether a synthetic data
generator can recover an unknown population distribution, but whether a
substantive conclusion remains stable under plausible perturbations of the
observed sample.

SYDNI addresses this narrower problem. The method generates synthetic
extensions of an observed dataset for use as secondary robustness checks. The
synthetic rows are deliberately not treated as new empirical evidence. Instead,
they are used to stress-test whether a downstream conclusion changes when the
observed data are resampled, perturbed, regularized, and passed through a
structure-preserving synthesis engine.

The `sydni` R package implements this workflow through a single front end,
`sydni()`, with five synthesis engines:

- `engine = "brm"`: conditional Bayesian regression synthesis using `brms`.
- `engine = "copula"`: subgroup-aware Gaussian copula synthesis with empirical
  margins.
- `engine = "bart"`: nonlinear conditional synthesis using Bayesian additive
  regression trees.
- `engine = "factor"`: subgroup-aware latent factor synthesis for psychometric
  item batteries.
- `engine = "graphical"`: subgroup-aware sparse Gaussian graphical synthesis.

The main contribution of SYDNI is not a claim that synthetic augmentation
increases evidential sample size. Its contribution is a conservative
perturbation framework for asking whether analytic conclusions are fragile to
sampling noise, modeling assumptions, subgroup imbalance, and dependence
structure. This positioning leads to three design principles:

1. **Structure preservation.** Synthetic rows should resemble the observed data
   in the aspects relevant to the planned analysis.
2. **Conservative perturbation.** Synthetic generation should not mechanically
   create stronger evidence by increasing row count.
3. **Robustness interpretation.** Augmented analyses should be reported as
   secondary sensitivity analyses, not as replacements for original-sample
   analyses.

The remainder of this paper draft describes the SYDNI workflow, the
engine-specific methods, the software interface, recommended validation, a
benchmarking design, and pilot-scale benchmark results from two simulation
studies - a manifest-variable inferential-safety benchmark (Section 8.6) and
a measurement-aware calibration benchmark (Section 8.7) - followed by a
single-run empirical illustration on a real psychometric dataset
(Section 9) that checks both studies' findings against real data.

## 2. Methodological premise

Let

\[
D = \{x_i\}_{i=1}^{n_{\mathrm{obs}}}
\]

denote the observed dataset. Let

\[
V = \{v_1, \ldots, v_p\}
\]

be the variables to synthesize and

\[
C = \{c_1, \ldots, c_q\}
\]

be optional covariates. Let \(n_{\mathrm{aug}}\) denote the requested synthetic
sample size and let \(D^\star\) denote the synthetic dataset returned by SYDNI.

SYDNI first constructs a feeder dataset from the observed data. The feeder data
are obtained by resampling observed rows and injecting controlled noise into
the variables to be synthesized. The noisy feeder data are then passed to a
synthesis engine. Each engine estimates a different representation of structure:
conditional regression structure, dense dependence structure, nonlinear
conditional structure, latent factor structure, or sparse direct dependence.

This design distinguishes SYDNI from two common interpretations of synthetic
data. First, SYDNI is not a privacy mechanism on its own. It may be possible to
evaluate disclosure risk for SYDNI outputs, but the package does not claim
formal privacy protection. Second, SYDNI is not an estimator that turns a small
sample into a larger effective sample. The synthetic data are generated from the
observed sample and should therefore be interpreted as sensitivity-analysis
material.

## 3. General workflow

At the top level, the user calls:

```r
result <- sydni(...)
```

The returned object has the form:

```r
list(
  models = ...,
  data = synthetic_data,
  meta = metadata
)
```

where `models` contains engine-specific fitted objects or summaries, `data`
contains the synthetic dataset, and `meta` stores lightweight generation
metadata. The original input dataset is not stored in `meta`.

### 3.1 Column resolution and validation

If `variables = NULL`, SYDNI synthesizes all non-covariate columns:

\[
V = \mathrm{names}(D) \setminus C.
\]

The implementation validates that `data` is a data frame or tibble, that
`variables` and `covariates` exist in `data`, and that the two sets do not
overlap. Synthesized columns must be numeric, integer, factor, or ordered
factor, and each synthesized column must have at least two distinct non-missing
values.

### 3.2 Synthetic sample size

The user specifies either `n_aug` directly or `prop_aug` proportionally. If
`prop_aug` is supplied, then:

\[
n_{\mathrm{aug}} =
\mathrm{round}(n_{\mathrm{obs}} \cdot \mathrm{prop\_aug}).
\]

If neither argument is supplied, the current implementation defaults
`n_aug` to the observed row count. If `replace = FALSE` but the requested
synthetic sample size exceeds the observed row count, the implementation sets
`replace = TRUE` with a warning.

### 3.3 General algorithm

```text
Algorithm 1: SYDNI general workflow

Input:
  Observed data D
  Synthesized variables V
  Optional covariates C
  Synthetic sample size n_aug
  Engine E
  Noise and shrinkage settings theta

1. Validate D, V, and C.
2. If n_aug is unspecified, set n_aug = n_obs.
3. Construct a feeder dataset by bootstrap resampling observed rows.
4. Inject controlled noise into synthesized variables.
5. Dispatch to engine E:
   a. If E is conditional, fit one model per target variable.
   b. If E is dependence-first, estimate subgroup dependence and margins.
   c. If E is latent-factor, estimate subgroup latent structure.
6. Generate synthetic draws from the fitted engine.
7. Reconstruct synthetic values to observed support and storage type.
8. Append or preserve covariates according to the engine design.
9. Return fitted objects, synthetic data, and lightweight metadata.
```

## 4. Noise-informed shrinkage

Noise injection is central to the SYDNI workflow. Its practical purpose is to
discourage trivial copying of observed rows and to weaken dependence on
sample-specific idiosyncrasies. Noise is applied by `inject_noise()` and
`add_noise()`.

A variable is treated as integer-like if it is stored as integer, factor,
ordered factor, or numeric with whole-number non-missing values. After noise is
added, values are clamped either to the observed range or to a user-supplied
`noise_bounds = c(lower, upper)`.

The supported noise modes are:

- `noise_method = "discrete"`: integer-like variables receive noise in
  \(\{-1, 0, +1\}\), with probabilities \(c/2\), \(1-c\), and \(c/2\);
  continuous variables receive uniform noise scaled by the observed standard
  deviation.
- `noise_method = "gaussian"`: variables receive Gaussian noise with scale
  \(c \cdot s_x\), where \(s_x\) is the observed standard deviation.
- `noise_method = "truncnormal"`: variables receive zero-centered truncated
  Gaussian noise bounded by \(\pm c \cdot s_x\).
- `noise_method = "random"`: variables receive uniform noise in `[-r, r]`,
  independent of observed scale.

When `c = "auto"`, discrete noise uses `p` as the probability control for a
non-zero perturbation, while continuous noise uses \(c = \sqrt{f}\), so that
the implied noise variance is approximately \(f\) times the observed variance.

Shrinkage is engine-specific. In the BRM engine, slope priors shrink as the
augmentation ratio grows:

\[
\sigma_\beta =
\frac{\mathrm{base\_sd}}
     {\left(1 + n_{\mathrm{aug}} / n_{\mathrm{obs}}\right)^k}.
\]

The BART engine increases tree regularization as \(n_{\mathrm{aug}} /
n_{\mathrm{obs}}\) grows. The copula engine shrinks subgroup correlations
toward a global correlation target. The factor engine shrinks loadings and
raises uniqueness floors. The graphical engine increases the effective EBIC
penalty for sparse precision selection. Across engines, the shared design goal
is to make larger synthetic augmentations more conservative, not more
evidentially aggressive.

### 4.1 Measurement-aware calibration

All five engines can additionally be run in a measurement-aware mode by
supplying a `measurement_blocks` argument, a `between_block_shrink` argument,
or both. This mode targets a narrower but common applied setting: item
batteries and composite scores, where the analyst wants within-scale
measurement structure protected from shrinkage (item margins, within-block
correlations, internal-consistency reliability, average variance extracted)
while between-construct relationships are made conservative as the synthetic
sample grows.

`measurement_blocks` accepts either a plain named list of item vectors, where
each entry is simultaneously a leaf measurement block and a final structural
family, or a hierarchical composite-list input in which lower-order item
batteries remain leaves and batteries referenced by a higher-order composite
are combined into one final structural family. Any synthesized variable not
covered by a supplied block becomes an implicit singleton family, so
supplying only `between_block_shrink` with no blocks activates uniform
shrinkage across every synthesized variable.

Measurement-aware calibration proceeds in two conceptual layers, implemented
differently across the five engines depending on whether the engine is a
correlation-matrix simulation engine (copula, graphical), a latent-factor
engine (factor), or a conditional-model engine (BRM, BART):

1. **Structure-level attenuation.** For the matrix-simulation and
   latent-factor engines, between-family dependence is attenuated by a
   multiplier \(\lambda \in [0,1]\) (`between_block_shrink`) while
   within-family structure is restored from the observed reference before
   simulation. The BRM engine instead fits models directly at a structural
   family-score layer, with between-family predictor counts and priors
   shrunk by the same multiplier. The BART engine has no structure-level
   attenuation step, since it is not a correlation-matrix-driven engine.

2. **Realized-score calibration.** After synthetic output is generated, let
   \(R_{\mathrm{obs}}\) be the observed correlation matrix of standardized
   structural-family scores. The calibration target is
   \[
   R_{\mathrm{target}} = I + \lambda (R_{\mathrm{obs}} - I).
   \]
   A deterministic Iman-Conover-style (Iman & Conover, 1982) rank-permutation search evaluates
   twenty-five candidate calibrations, each applying one row permutation to
   every item within a structural family so that each family's complete
   multivariate row distribution moves as a unit. The candidate with the
   smallest maximum off-diagonal target error is selected, using mean
   absolute error as a tie-breaker. Because the uncalibrated output is
   always included as a candidate, calibration cannot produce worse target
   error than the raw generated data, and because complete families move
   together, within-family correlations, item margins and support, and
   reliability statistics such as alpha, composite reliability, and average
   variance extracted are unaffected by calibration.

The BRM engine applies realized-score calibration before item-level
reconstruction rather than after, because its item values are reconstructed
by copying complete observed item blocks from donor rows matched on the
calibrated family score. This preserves valid response support exactly but
means BRM's final item-level composite correlations can differ slightly from
the calibration target, unlike the other four engines, which calibrate
already-generated item-level data directly.

### 4.2 Retention-aware augmentation and group balancing

SYDNI also supports generating a larger internal candidate pool and
returning only a retained subset matching a target covariate value, through
two mutually exclusive interfaces: `retain_by`/`retain_value`, and
`balance_by`/`balance_value`/`balance_target = "largest_group"`. The
balancing interface computes the required augmentation automatically as the
difference between the largest observed group size and the observed size of
the target group, then delegates to the same retention mechanism.

Given a requested retained count \(n_{\mathrm{retained}}\) and the observed
prevalence \(\pi\) of the retained group in the original data, the internal
candidate pool size is
\[
n_{\mathrm{candidate}} =
\begin{cases}
n_{\mathrm{retained}}, & \pi = 1,\\[4pt]
n_{\mathrm{retained}} +
Q_{\mathrm{NB}}(0.999;\ n_{\mathrm{retained}},\ \pi), & \pi < 1,
\end{cases}
\]
where \(Q_{\mathrm{NB}}\) is the negative-binomial quantile function. This
sizes the candidate pool so that, with 99.9% confidence, at least
\(n_{\mathrm{retained}}\) generated rows will match the target group.

Critically, the chosen engine receives the inflated candidate-pool size as
its generation target but receives the originally requested retained count
separately for augmentation-linked shrinkage and regularization. This
decouples shrinkage strength from the size of the internal candidate pool: a
larger candidate pool built to reliably retain a small target group does not
make the returned synthetic rows more or less conservative than requesting
that same retained count directly. After generation, eligible rows are
identified, a deterministic sample of the requested size is drawn without
replacement, and the retained subset is passed through the same
realized-score calibration described in Section 4.1, targeting a correlation
derived from the observed structural-score correlation in the full original
dataset. This recalibration step runs even when the chosen engine has no
native measurement-aware path.

When `between_block_shrink = "auto"` (or omitted while measurement blocks or
retention are in use), shrinkage is computed as
\[
\lambda_{\mathrm{auto}} =
\frac{1}{1 + n_{\mathrm{retained}} / n_{\mathrm{obs}}},
\]
based on the number of synthetic rows returned for analysis rather than the
size of any internal candidate pool.

## 5. Engine architecture

The five engines are best understood as three methodological families.

### 5.1 Conditional engines

The conditional engines synthesize each target variable from selected
predictors and covariates.

The BRM engine fits one `brms` model per synthesized variable. Predictor sets
can be full, random, information-ranked, or randomly sampled from
information-ranked candidates. Families can be supplied manually or inferred
from the variable type. The fitted models generate posterior predictive draws
for feeder rows, and outputs are reconstructed to the original variable support.

The BART engine follows the same broad target-wise structure but replaces
parametric Bayesian regressions with Bayesian additive regression trees through
`dbarts`. This engine is intended for nonlinear conditional structure and
interactions. Binary targets are reconstructed by thresholding; ordinal-like
targets are handled through a numeric approximation path and then mapped back
to observed support.

### 5.2 Dependence-first engines

The dependence-first engines model joint structure rather than fitting one
conditional model per target.

The copula engine defines subgroups from covariate combinations, constructs
noisy subgroup feeder data, estimates subgroup Gaussian copula correlations,
shrinks subgroup correlations toward a global correlation target, and maps
simulated Gaussian copula uniforms through empirical subgroup margins.

The graphical engine also works subgroup-wise, but it targets sparse direct
dependence. It transforms subgroup feeder data into a latent Gaussian space,
fits a sparse precision matrix through graphical lasso using `huge`, converts
the selected precision matrix to a stabilized correlation matrix, simulates
Gaussian draws, and reconstructs observed margins empirically. If sparse
estimation fails, it falls back to a stabilized latent correlation matrix.

### 5.3 Latent-structure engine

The factor engine is intended for ordinal or integer-like item batteries. It
constructs noisy subgroup feeder data, identifies non-constant active items,
estimates a correlation matrix, fits a factor model using `psych`, shrinks
loadings and factor correlations, simulates latent item values, and maps those
values back through empirical subgroup margins. If the factor model is not
viable in a subgroup, the implementation falls back to empirical subgroup
reconstruction.

### 5.4 Engine comparison

| Engine | Statistical target | Covariate handling | Main strength | Main caution |
|---|---|---|---|---|
| `brm` | Conditional regression | Native predictors | Explicit modeled relationships and mixed outcome support; never the worst engine on any manifest-variable safety axis (Section 8.6) | Computational cost and model specification; in measurement-aware mode, pilot benchmarking (Section 8.7) found no quality advantage over the cheaper engines large enough to offset a further 400-1,400x runtime premium |
| `bart` | Nonlinear conditional regression | Native predictors | Nonparametric flexibility; theoretically suited to nonlinearities and interactions; on manifest variables, a near-free `brm` stand-in (Section 8.6) | Ordinal targets use an approximation path; a direct oracle-model test (Section 8.6) found no measurable nonlinear/interaction advantage over the other four engines; in measurement-aware mode (Section 8.7) its within-scale item-reliability retention is 3-4x worse than the other engines, though it leads on power recovery and between-construct fidelity there |
| `copula` | Dense joint dependence | Subgrouping | Empirical margins and dependence preservation; fast default alternative on manifest variables (Section 8.6) | Tail behavior may be smoothed |
| `factor` | Latent measurement structure | Subgrouping | Psychometric item batteries; best-balanced measurement-aware engine in Section 8.7's pilot | Sparse subgroups can be fragile; lowest power recovery among measurement-aware engines in Section 8.7's pilot |
| `graphical` | Sparse direct dependence | Subgrouping | Network-like conditional structure; fastest and near-best-balanced measurement-aware engine in Section 8.7's pilot | Sparsity must be scientifically meaningful |

All five engines additionally accept the measurement-aware calibration and
retention-aware augmentation controls described in Sections 4.1 and 4.2. The
underlying algorithmic target in the table above (conditional regression,
dense dependence, and so on) does not change; measurement-aware mode adds a
between-construct shrinkage and calibration layer on top of each engine's
existing generative mechanism. Engine recommendations are scenario-dependent
rather than a single global ranking - see Section 8.7's summary table for
how the recommended engine differs depending on whether the analyst has
manifest/composite variables only, has indicator-level data but only needs
composite-score output, or needs the synthesized indicators themselves.
Full step-by-step pseudocode for all five engines is given in Appendix A.

## 6. Software interface

The public package interface centers on `sydni()`:

```r
sydni(
  data,
  variables = NULL,
  covariates = NULL,
  engine = c("brm", "copula", "bart", "factor", "graphical"),
  noise_method = c("truncnormal", "gaussian", "discrete", "random"),
  c = "auto",
  round_int = TRUE,
  p = 0.5,
  f = 0.5,
  r = 1,
  noise_bounds = NULL,
  predictor_method = c("random_information", "random", "information", "full"),
  predictor_k = NULL,
  replace = TRUE,
  n_aug = NULL,
  prop_aug = NULL,
  retain_by = NULL,
  retain_value = NULL,
  balance_by = NULL,
  balance_value = NULL,
  balance_target = "largest_group",
  families = NULL,
  base_sd = 10,
  sd_scale = 2,
  k = 2,
  chains = 2,
  cores = "auto",
  iter = 500,
  warmup = floor(iter / 2),
  thin = 2,
  algorithm = c("fullrank", "meanfield", "sampling", "pathfinder", "laplace", "fixed_param"),
  backend = c("rstan", "cmdstanr"),
  engine_control = NULL,
  measurement_blocks = NULL,
  between_block_shrink = NULL,
  seed = NA
)
```

`measurement_blocks` and `between_block_shrink` activate the measurement-aware
calibration described in Section 4.1 and are supported by all five engines.
`retain_by`/`retain_value` and `balance_by`/`balance_value`/`balance_target`
activate the retention-aware augmentation and group-balancing behavior
described in Section 4.2 and are mutually exclusive with each other and,
in the case of `balance_by`, with directly supplying `n_aug` or `prop_aug`.

The engine functions are also exported:

- `engine_brm()`
- `engine_copula()`
- `engine_bart()`
- `engine_factor()`
- `engine_graphical()`

Noise utilities are exported as `inject_noise()` and `add_noise()`.

### 6.1 Example using the included grit data

The package includes `grit`, a dataset with 500 OpenPsychoMetrics responses to
Big Five personality items and the Grit Scale. The following example generates
synthetic data for openness and extraversion items while preserving two grit
items as covariates:

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

```
  o1 o2 o3 o4 o5 o6 o7 o8 o9 o10 e1 e2 e3 e4 e5 e6 e7 e8 e9 e10 gs1 gs2
1  2  5  3  3  2  5  4  2  4   5  3  2  2  3  5  2  3  2  3   1   5   2
2  5  4  4  5  2  5  3  4  4   4  2  5  3  4  4  4  4  3  4   3   5   2
3  3  4  5  5  4  5  3  3  5   5  1  3  4  2  3  3  4  2  3   1   5   2
4  2  3  4  4  4  3  4  3  4   5  3  3  3  3  4  5  1  1  1   4   5   2
5  2  2  4  2  3  2  4  3  4   3  1  4  3  1  4  1  5  1  3   1   5   2
6  4  5  5  5  1  5  4  5  4   5  1  1  4  1  4  3  3  1  4   2   5   2
```

`result$data` has the same columns and support as `grit_complete` (Likert
items on their original 1-5 scale, plus the two preserved covariates).
`gs1` and `gs2` retain real observed Grit-item values throughout, since
covariates are drawn from bootstrap-resampled observed rows rather than
being synthesized or noise-injected (Section 3).

The same workflow can be run with another engine by changing only the `engine`
argument and, when needed, engine-specific controls:

```r
engines <- c("copula", "factor", "graphical")

outputs <- lapply(engines, function(e) {
  sydni(
    data = grit_complete,
    variables = c(sprintf("o%01d", seq(10)), sprintf("e%01d", seq(10))),
    covariates = c("gs1", "gs2"),
    n_aug = 100,
    engine = e,
    seed = 123
  )
})

names(outputs) <- engines
```

Conditional engines can be used through the same interface:

```r
bart_out <- sydni(
  data = grit_complete,
  variables = c(sprintf("o%01d", seq(10)), sprintf("e%01d", seq(10))),
  covariates = c("gs1", "gs2"),
  n_aug = 100,
  engine = "bart",
  predictor_k = 5,
  engine_control = list(num_trees = 50, num_burn = 250, num_draws = 250),
  seed = 123
)
```

Measurement-aware calibration can be requested for any engine by declaring
item blocks and a shrinkage strength. The following call protects the
openness and extraversion item batteries as separate structural families
while shrinking their between-family correlation toward half its observed
value:

```r
measurement_out <- sydni(
  data = grit_complete,
  variables = c(sprintf("o%01d", seq(10)), sprintf("e%01d", seq(10))),
  covariates = c("gs1", "gs2"),
  n_aug = 100,
  engine = "copula",
  measurement_blocks = list(
    openness = sprintf("o%01d", seq(10)),
    extraversion = sprintf("e%01d", seq(10))
  ),
  between_block_shrink = 0.5,
  seed = 123
)
```

Comparing observed and synthetic correlation structure directly shows the
two-layer mechanism from Section 4.1 at work on this single call: the
within-block (openness item) correlation shrinks only modestly, while the
between-block (openness-extraversion composite) correlation shrinks far
more sharply, roughly in line with the requested 0.5 shrinkage factor
applied to a modest observed association. This single-call comparison is
illustrative rather than a validated benchmark result; Section 8.7
reports the same flat-versus-measurement-aware contrast at full benchmark
scale.

```r
o_items <- sprintf("o%01d", seq(10))
e_items <- sprintf("e%01d", seq(10))

obs_within <- mean(cor(grit_complete[o_items])[lower.tri(diag(10))])
syn_within <- mean(cor(measurement_out$data[o_items])[lower.tri(diag(10))])

obs_between <- cor(rowMeans(grit_complete[o_items]), rowMeans(grit_complete[e_items]))
syn_between <- cor(rowMeans(measurement_out$data[o_items]), rowMeans(measurement_out$data[e_items]))
```

```
Mean within-openness item correlation: observed = 0.284, synthetic = 0.183
Openness-extraversion composite correlation: observed = 0.174, synthetic = 0.010
```

Retention-aware augmentation instead targets a specific covariate-defined
group. For an observed sample containing an imbalanced grouping covariate,
the following call returns enough synthetic rows of the underrepresented
group to match the largest observed group, with automatic measurement-aware
shrinkage applied to the retained subset:

```r
balanced_out <- sydni(
  data = dat,
  variables = item_vars,
  covariates = "group",
  engine = "copula",
  balance_by = "group",
  balance_value = "minority_level",
  between_block_shrink = "auto",
  seed = 123
)
```

## 7. Validation and reporting

SYDNI output should be validated before being used in a downstream robustness
analysis. At minimum, analysts should compare:

- marginal distributions of synthesized variables;
- subgroup proportions induced by covariates;
- correlation, rank-correlation, or other dependence summaries;
- downstream model coefficients, predictions, classifications, or decisions;
- sensitivity to `engine`, `n_aug`, `noise_method`, `p`, `f`, `c`, `r`, and
  `k`.

The original-sample analysis should be reported first. SYDNI analyses should be
reported alongside the original analysis, not in place of it. If a result is
present in both the original and augmented analyses, the agreement may be
reassuring. If a result is absent in the original data but appears only after
augmentation, it should be interpreted as fragile, suggestive, or exploratory
rather than confirmatory.

Recommended reporting language:

> SYDNI analyses were treated as secondary robustness checks alongside the
> original-sample analyses. Synthetic augmentation was not interpreted as
> additional empirical evidence, and effects that appeared only after
> augmentation were interpreted cautiously.

Manuscripts using SYDNI should report the full generation settings, including
`engine`, `n_aug` or `prop_aug`, `noise_method`, noise controls, predictor
selection settings, shrinkage settings, backend settings, and seed.

## 8. Benchmarking design

The benchmark design evaluates the claim SYDNI is actually built to support:
conservative inferential safety under synthetic data augmentation. The goal
is not to show that SYDNI is superior to SMOTE or ADASYN for classifier
balancing, superior to GANs for synthetic training data, or superior to
`synthpop` for disclosure-oriented synthetic microdata. Instead, these
methods serve as off-label comparators, included because analysts may be
tempted to reach for them for sample augmentation in social and behavioral
analyses even though that is not their intended use.

The primary analysis target is the pooled augmented dataset:

\[
D_{\mathrm{pooled}} = D_{\mathrm{obs}} \cup D_{\mathrm{syn}}.
\]

Synthetic-only analyses are retained as diagnostics, but the core question
each benchmark asks is whether adding synthetic rows makes ordinary
statistical analyses look stronger, more precise, or more significant than
the observed data justify.

### 8.1 Simulation regimes

The benchmark is simulation-first because only simulations provide known
truth. Both pilots (Sections 8.6-8.7) used:

\[
n_{\mathrm{obs}} \in \{75, 150, 300\},
\qquad
n_{\mathrm{syn}} / n_{\mathrm{obs}} \in \{0.5, 1, 2\}.
\]

Replication count was chosen relative to what each comparison needed to
resolve rather than fixed to an arbitrarily large target: a method-family
gap large enough to matter substantively is already resolved by Monte
Carlo error at a modest replication count, while a close within-family
comparison is not made trustworthy by more replications of the same
design, since its dominant uncertainty is external validity - how well
the simulated design matches real data - not sampling noise. Sections
8.6-8.7 accordingly report pilot results at 10 replications per condition,
treating large directional gaps as resolved and close comparisons as
directional rather than confirmed; see the discussion in Section 12.

Section 8.6's manifest-variable pilot covers seven simulation scenarios,
chosen to isolate distinct failure modes:

1. **Null-only inference.** Many predictors and outcomes with no true effects.
   This is the cleanest Type I error test.
2. **Sparse true effects.** A small number of true predictors among many null
   predictors, testing false discovery and effect inflation.
3. **Correlated predictors.** Latent-factor correlated predictors with true and
   null variables correlated, testing spuriousness under confounding-like
   structure.
4. **Ordinal psychometric data.** Likert items generated from latent factors,
   testing polychoric/correlation recovery, factor loading recovery, and
   reliability recovery.
5. **Moderation and subgroup imbalance.** Imbalanced group variables with true
   and null interaction terms, testing augmentation-induced moderation.
6. **Nonlinear structure.** Thresholds or interactions, giving flexible methods
   such as BART a fair setting.
7. **Small-cell stress test.** Covariate-defined cells with very small observed
   counts, testing subgroup instability.

Section 8.7's companion pilot replaces these seven scenarios with a single
latent-variable item-battery design purpose-built to exercise the
measurement-aware mechanism (Section 4.1), described in full there.

### 8.2 Methods compared

Each pilot compared:

- observed-only analysis as the reference baseline;
- bootstrap-only and noise-only bootstrap as minimal baselines;
- SMOTE and ADASYN from `smotefamily`, used off-label with a binary balancing
  target;
- `synthpop`, in both a sequential/CART configuration and a conservative
  sample-based configuration;
- all five SYDNI engines (copula, graphical, BART, factor, BRM) - full BRM was
  included in both pilots at its real computational cost rather than a
  reduced subset, since runtime, while substantial, remained feasible.

A tabular GAN comparator (e.g. CTGAN) was considered but not included in
either pilot, since a reproducible, environment-independent comparison could
not be guaranteed; Section 8.6 notes this exclusion explicitly. Failures,
warnings, degenerate outputs, constant columns, and dropped model terms are
reported as benchmark outcomes rather than silently discarded.

### 8.3 Primary inferential-safety metrics

The headline metrics directly test whether augmentation distorts inference:

- Type I error rate for true-null coefficients.
- False discovery rate across many null predictors.
- Augmentation-only significance rate: nonsignificant in observed data but
  significant after augmentation.
- Effect-size inflation relative to truth:
  \[
  |\hat{\beta}_{\mathrm{aug}}| - |\beta_{\mathrm{true}}|.
  \]
- Effect-size inflation relative to observed-only estimates:
  \[
  |\hat{\beta}_{\mathrm{aug}}| - |\hat{\beta}_{\mathrm{obs}}|.
  \]
- Standard-error compression:
  \[
  \mathrm{SE}_{\mathrm{aug}} / \mathrm{SE}_{\mathrm{obs}}.
  \]
- 95% interval coverage against known truth.
- Sign reversal rate for true nonzero effects.
- Moderation false-positive rate for true-null interactions.

### 8.4 Secondary structure and authenticity metrics

Structure-preservation metrics explain why inferential distortion may occur:

- Mean absolute correlation matrix error.
- Polychoric correlation error for Likert items.
- False positive and false negative correlation-edge rates.
- Factor loading RMSE and loading sign agreement.
- Alpha/omega recovery error.
- Likert category distribution error.
- Subgroup proportion error.

Diagnostic authenticity metrics include duplicate and near-duplicate rates,
support violations, impossible category rates, constant-column rates,
dropped-term rates, warning counts, failure rates, runtime, and memory cost.

### 8.5 Result displays

This paper does not present a single overall winner; results are instead
ranked by benchmark family, using the following displays:

- **Table 1:** Type I error, false discovery, coverage, and augmentation-only
  significance by scenario and method.
- **Table 2:** Effect inflation and standard-error compression by augmentation
  ratio.
- **Table 3:** Correlation, psychometric, and subgroup-structure preservation
  metrics.
- **Figure 1:** Type I error, false discovery rate, and coverage by method
  across sample sizes.
- **Figure 2:** Augmentation-only significance rate and power recovery rate
  across sample sizes.
- **Figure 3:** Runtime diagnostics by method.
- **Figure 4:** Flat-vs-measurement-aware paired comparison by engine across
  sample sizes.

Section 8.6 reports pilot-scale numbers for Tables 1-3 and Figures 1-3 on
manifest variables; Section 8.7 reports a companion pilot isolating the
measurement-aware calibration mechanism on a latent-variable item-battery
design, including Figure 4; Section 9 reports a single-run empirical
check of both against real data. See Section 12 for why additional
replications of the existing designs are not treated as a priority
alongside the results already reported here.

### 8.6 Pilot results

A pilot run of the design above has been completed and is reported here.
It covers all 7 scenarios from Section 8.1, the
full \(n_{\mathrm{obs}} \in \{75, 150, 300\}\) and augmentation-ratio grid
from Section 8.1, 10 replications per condition, 12 methods spanning every
category in Section 8.2 except a tabular GAN comparator, and 7,560 total
simulated conditions. Full results, per-scenario detail, and reproducing
scripts are reported in `results/benchmarks/FINDINGS.md` in the package
repository; this section summarizes the headline findings.

**Off-label stress test.** Every SYDNI engine and `synthpop`'s CART
variant clearly outperformed the off-label bootstrap/SMOTE/ADASYN
comparators on Type I error and false discovery rate, by roughly a factor
of two, consistently across all three sample sizes tested. This was true
even though these comparator methods were evaluated under the same
discouraged pooled-analysis usage as SYDNI itself in this stress test, not
under their own intended use cases.

![Figure 1: Off-label stress test - Type I error rate, false discovery rate, and 95% coverage across n_obs, by method.](figures/figure1_type1_coverage.png)

**Figure 1.** Off-label stress test: Type I error rate, false discovery
rate, and 95% coverage across \(n_{\mathrm{obs}}\), for the five SYDNI
engines (left) against the five off-label comparators evaluated under the
same discouraged pooled-analysis usage (right). The gray band and dashed
line in each panel trace the min-max range and mean of the *other* group
at each \(n_{\mathrm{obs}}\), so a SYDNI engine's line can be read
directly against the comparator spread (and vice versa) without
cross-referencing the opposite panel. Method labels sit to the left of
\(n_{\mathrm{obs}} = 75\) at each line's true starting value, nudged
vertically only when two methods would otherwise overlap. `observed_only`
and `synthpop`'s sample-based variant are omitted from the figure - the
former is not an augmentation method and the latter is a documented floor
baseline rather than a real comparator (see the pilot results below) -
but both remain in the underlying pooled results summarized in this
section.

**Intended-use safety.** Under SYDNI's documented secondary-check
workflow, the false-positive rate for effects that were null in the
observed data but spuriously significant only after augmentation
(`false_augmentation_only_sig_rate`) was lowest for `sydni_factor`
(0.089) and `sydni_graphical` (0.096), and highest among the SYDNI/synthpop
engines for `sydni_bart` (0.124) - all still well below every off-label
comparator (0.151-0.185).

![Figure 2: Intended-use safety - false augmentation-only significance rate and power recovery rate across n_obs, by method.](figures/figure2_aug_only_sig.png)

**Figure 2.** Intended-use safety: false augmentation-only significance
rate (top) and power recovery rate (bottom) across
\(n_{\mathrm{obs}}\), same panel/label/reference-band conventions as
Figure 1. These are deliberately not the same quantity as Figure 1's
false discovery rate: `false_augmentation_only_sig_rate` is restricted to
terms with no true effect that were specifically non-significant in the
observed-only fit and became significant only after augmentation - the
"manufactured finding" SYDNI's documented secondary-check workflow is
meant to catch - while `power_recovery_rate` is the same
observed-to-augmented transition restricted to terms with a real, nonzero
effect, a positive event rather than a safety risk. False discovery rate
(Figure 1) is broader still: the share of *all* significant results in
the augmented fit that are false, regardless of whether they were already
significant before augmentation.

**Structure retention and engine tradeoffs.** No single SYDNI engine
dominated on every axis; the five engines were designed as different
tradeoffs rather than competitors on one scale, and the pilot's role is to
show what each tradeoff actually costs and buys.

\begin{center}
\begin{tabular}{lrrrr}
\hline
Engine & Type I error (rank) & Power recovery (rank) & Correlation error (rank) & sec/condition \\
\hline
\texttt{sydni\_factor}    & 0.135 (1st, safest) & 0.084 (3rd)          & 0.082 (5th, worst)  & 1.88 \\
\texttt{sydni\_graphical} & 0.145 (2nd)         & 0.077 (5th, worst)   & 0.077 (4th)         & 0.04 \\
\texttt{sydni\_brm}       & 0.151 (3rd)         & 0.082 (4th)          & 0.056 (3rd)         & 96.8 \\
\texttt{sydni\_bart}      & 0.173 (4th)         & 0.085 (2nd)          & 0.053 (2nd)         & 0.39 \\
\texttt{sydni\_copula}    & 0.178 (5th, worst)  & 0.099 (1st, best)    & 0.039 (1st, best)   & 0.04 \\
\hline
\end{tabular}
\end{center}

**`sydni_factor`'s row is not directly comparable in magnitude to the
other four engines' rows.** `engine_factor()` requires ordinal/integer
input, so its numbers come from a discretized version of the continuous
outcome `y`, while the other four engines' numbers come from the same
continuous `y` used elsewhere in this pilot. `sydni_factor`'s rank
positions above are directional evidence about the factor engine on its
own terms, not a like-for-like comparison against the other four -
"1st, safest" and similar labels should be read with that qualification
rather than as a claim that `sydni_factor` is safer than the other
engines under an identical outcome.

Individual rank order among these five engines should be read as
directional at this pilot's 10 replications per condition rather than as
a confirmed ranking (Section 8.6 continues below); the large gap between
the SYDNI/synthpop family and the off-label comparators above is far more
robust to replication count than the closer distinctions within this
table. `sydni_brm` is not the power-recovery leader among SYDNI engines,
contrary to an earlier draft of this analysis based on the
\(n_{\mathrm{obs}} \in \{75, 150\}\)-only pilot: that ranking flipped once
\(n_{\mathrm{obs}} = 300\) was added, and \texttt{sydni\_copula} leads on
power recovery in the full pooled data.

The runtime column changes the interpretation of this table substantially.
`sydni_brm` is the one engine that is never worst on any axis (3rd, 4th,
2nd) - a real, defensible sense of "balanced" - but costs roughly
250-2,500x more compute than every other engine. Comparing each engine
against `sydni_brm` directly (gap = engine minus `sydni_brm`; negative
Type I error/correlation gap and positive power gap both indicate the
cheaper engine performing better) shows how much of that balance survives
at a fraction of the cost:

\begin{center}
\begin{tabular}{lrrrr}
\hline
Engine & Speedup vs.\ \texttt{sydni\_brm} & Type I gap & Power gap & Correlation gap \\
\hline
\texttt{sydni\_graphical} & 2{,}464x & $-$0.006 (safer)  & $-$0.004 (\textasciitilde tied) & $+$0.021 (worse) \\
\texttt{sydni\_copula}    & 2{,}166x & $+$0.027 (worse)  & $+$0.017 (better)               & $-$0.017 (better) \\
\texttt{sydni\_bart}      & 250x     & $+$0.022 (worse)  & $+$0.003 (\textasciitilde tied)  & $-$0.003 (\textasciitilde tied) \\
\texttt{sydni\_factor}    & 52x      & $-$0.016 (safer)  & $+$0.002 (\textasciitilde tied)  & $+$0.026 (worse) \\
\hline
\end{tabular}
\end{center}

This splits the family into three roles rather than a single ranking.
`sydni_brm` is the reference point when compute is not a binding
constraint and no single axis dominates the analyst's priority, at a cost
of the rare small/sparse-data ADVI fitting failure noted below.
`sydni_copula` is the most likely default fast alternative for most cases:
it posts the best power recovery and best structure retention of any
engine at a small fraction of `sydni_brm`'s compute cost, at the price of
the largest Type I error gap among the cheap alternatives. `sydni_bart`
and `sydni_graphical` are secondary, structurally different cheap
alternatives worth retaining rather than a first recommendation - each
also matches `sydni_brm` closely on two of three axes at 250-2,464x less
compute, giving up ground on a different single axis each. Neither engine
needs to lead any metric to be useful: no simulation-based benchmark can
realistically cover the full space of real-world data an analyst might
encounter, and having several structurally different, inexpensive options
available is itself part of the value of SYDNI's multi-engine design,
independent of any one engine's ranking in this pilot. `sydni_factor` is
the opposite-end specialist from `sydni_copula` rather than a `sydni_brm`
approximation: best safety in the family at the cost of the worst
structure retention, with power roughly neutral - the engine to reach for
when false-positive control is the explicit priority, subject to the
discretized-outcome caveat noted directly above the engine-tradeoff table.

Separately from the cross-engine comparison above, `sydni_brm` showed a
rare ADVI fitting failure confirmed confined to the smallest-\(n_{\mathrm{obs}}\),
largest-augmentation-ratio corner of the grid (4 failures out of 630
`sydni_brm` conditions, all at \(n_{\mathrm{obs}} = 75\), \(n_{\mathrm{aug}}/
n_{\mathrm{obs}} = 2\); 2 of the 4 in the small-cell scenario, 1 each in the
null-only and moderation/imbalance scenarios - see `KNOWN_ISSUES.md` in the
package repository), at a computational cost roughly two orders of
magnitude larger than any other engine (approximately 90 seconds per
condition versus fractions of a second).

![Figure 3: Mean runtime per condition by method, log scale.](figures/figure3_diagnostics.png)

**Figure 3.** Mean runtime per condition by method (log10 scale),
SYDNI engines against off-label comparators. `observed_only` and
`synthpop`'s sample-based variant are omitted as near-zero-runtime
reference/floor baselines rather than real comparators (see above). 4 of
7,560 conditions failed to complete (all `sydni_brm`, \(n_{\mathrm{obs}} =
75\); see `KNOWN_ISSUES.md` in the package repository). Every method also
emits occasional informational package notices (e.g., ordinal-as-numeric,
small-subgroup, or missing-`measurement_blocks` notices) unrelated to
fitting problems; the raw counts are reported alongside the reproducing
scripts in `results/benchmarks/FINDINGS.md`.

**Sample-size dependence.** Two patterns emerged from varying
\(n_{\mathrm{obs}}\) that were not previously anticipated in this design.
First, the power-recovery rate declined toward zero for every method as
\(n_{\mathrm{obs}}\) increased, consistent with augmentation mattering
primarily when the observed sample is itself underpowered - a pattern
that supports rather than undermines SYDNI's intended positioning.
Second, Type I error rate for null coefficients increased somewhat with
\(n_{\mathrm{obs}}\) for every method, including the observed-only
baseline that never touches augmentation; this likely reflects
main-effects-only analysis-model misspecification against a subset of
deliberately nonlinear/interactive simulation scenarios rather than an
augmentation-specific problem, and is a caveat on the absolute Type I
error figures rather than on the relative method rankings, which were
stable across sample sizes.

**Does BART show a distinct nonlinear/interaction advantage? A direct
test found none.** `sydni_bart`'s practical role - a structurally
different, near-free alternative to `sydni_brm` - is addressed above. A
narrower question remained open through two earlier follow-ups (relaxing
the engine's internal regularization settings, and a purpose-built binary
threshold-and-interaction scenario outside the main design), neither of
which surfaced BART's theoretically expected advantage: nonparametric
conditional-mean estimation should handle nonlinearity and interactions
without pre-specifying them (Section 5), in principle producing more than
parity with `sydni_brm`, not simply a cheaper match to it. Both follow-ups
left one hypothesis untested - that the benchmark's linear-main-effects
analysis model (\(y \sim x_1 + \cdots + x_{12} + \mathrm{group}\)) is
structurally incapable of rewarding correct nonlinear modeling regardless
of how faithfully an engine's synthetic data reproduces it.

A direct test of that hypothesis regenerated all 1,080 conditions from the
`nonlinear` scenario (\(y = 0.60 \cdot \mathbb{1}(x_1 \ge 4) + 0.45 \cdot
x_2 \cdot \mathbb{1}(x_3 \ge 4) + \varepsilon\)) from identical seeds and
re-analyzed them with an oracle-form model that includes the true terms
explicitly - \(y \sim \mathbb{1}(x_1 \ge 4) + x_2 \times \mathbb{1}(x_3
\ge 4) + x_4 + \cdots + x_{12} + \mathrm{group}\) - so an engine that
faithfully reproduced the nonlinear structure would receive full credit
for it. `sydni_bart` still showed no advantage: combined recovery of the
threshold and interaction terms placed `sydni_copula` (0.211) and
`sydni_bart`/`sydni_factor` (0.200/0.200) in an effective tie for best of
the five engines, with `sydni_bart`'s estimation error on both terms the
worst or second-worst of the five. This closes off the earlier
hypothesis rather than confirming it - an analysis model built
specifically to reward correct nonlinear/interaction recovery still finds
nothing distinguishing `sydni_bart` from the other four engines, so "the
metric could not see it" does not explain the earlier null results.
Two caveats limit how far this should be read: the naive off-label
comparators scored higher still on this metric, but only via the same
pseudo-replication artifact documented throughout this benchmark
(duplicating real rows trivially preserves the true pattern while also
inflating Type I error), not genuine structure recovery; and the true
effects were already significant in the observed-only fit 60-64% of the
time before augmentation, making the recovery event rare enough that the
narrow 0.156-0.211 spread among the five SYDNI engines should not itself
be read as a confirmed ranking. The revised conclusion is the stronger,
more useful negative: no evidence of a `sydni_bart` nonlinear/interaction
advantage was found under a metric specifically designed to detect one,
plausibly because a threshold-and-interaction effect on a
Likert-discretized predictor still carries substantial monotonic signal
that correlation-preserving engines can partially capture without any
nonlinear machinery at all.

### 8.7 Measurement-aware calibration pilot results

Section 8.6's pilot tests inference on manifest predictor variables
directly - the variables synthesized are the variables analyzed. It never
exercises the measurement-aware calibration mechanism of Section 4.1
(`measurement_blocks` / `between_block_shrink`), which is specifically
designed to protect within-scale item reliability while shrinking
between-construct associations. A second, companion pilot was run to
isolate that mechanism's contribution directly, using a latent-variable
item-battery design: 4 predictor latent constructs (3, 5, 4, and 4
indicators respectively) and 1 outcome latent construct (4 indicators),
16 + 4 manifest items total, with the true structural effect placed on 2
of the 4 predictor constructs and randomized per replicate which two carry
a nonzero effect - decoupling which predictors a fixed block-list order
would select from which predictors are actually true.

Each of the five SYDNI engines was run twice per condition as a paired
ablation: once with no `measurement_blocks` argument (`flat`), and once
with the full measurement-aware path (`ma`). Critically, the downstream
analysis is identical either way - item-level data is always collapsed to
construct-level composite scores and analyzed as
\(Y_{\mathrm{score}} \sim F1_{\mathrm{score}} + F2_{\mathrm{score}} +
F3_{\mathrm{score}} + F4_{\mathrm{score}}\) - so the `flat` vs. `ma`
comparison isolates what the calibration mechanism itself contributes
during synthesis, independent of how the result is analyzed afterward.
The same 7 off-label/baseline comparators from Section 8.2 were included
for reference. Using the same \(n_{\mathrm{obs}} \in \{75, 150, 300\}\)
and augmentation-ratio grid as Section 8.1 at 10 replications per
condition, this produced 1,530 total conditions (17 methods: 7
comparators plus 5 engines \(\times\) 2 variants) with zero failures.
Full results and reproducing scripts are reported in
`results/benchmarks/simulation_study_2/FINDINGS.md` in the package
repository; this section summarizes the headline findings, organized
around three scenarios this design speaks to.

**Scenario 1: manifest/composite predictors with no indicator-level data.**
Not tested by this pilot - this is Section 8.6's scenario, and nothing
here revises its conclusions. The DGP and analysis target differ enough
(no latent/indicator/aggregation layer at all) that this pilot's numbers
do not generalize to that case in either direction.

**Scenario 2: indicator-level data exists, but only composite-score output
is needed.** Here the pilot's own flat-vs-ma comparison is directly
informative, and the result holds across every engine, not just one: every
`flat` variant's Type I error and false discovery rate (0.14-0.19 /
0.08-0.12) sits in the same range as the off-label comparators, while
every `ma` variant (0.03-0.05 / 0.01-0.03) sits close to the observed-only
baseline (0.028 / 0.020). The practical implication is that
`measurement_blocks` is worth supplying whenever indicator-level data
exists, even when the end goal is only a composite-score regression -
skipping it costs real inferential safety on every engine tested, at a
runtime premium that stays under a second per condition for
`sydni_copula`, `sydni_graphical`, and `sydni_bart`, but rises to roughly
20x for `sydni_brm` (28.8 to 585.9 seconds).

![Figure 4: Flat vs. measurement-aware, paired by engine, across n_obs.](figures/figure4_flat_vs_ma.png)

**Figure 4.** Type I error rate, false discovery rate, and power recovery
rate for each of the five SYDNI engines run without `measurement_blocks`
(`flat`, left) against the same five engines run with the full
measurement-aware calibration path (`ma`, right), paired by engine and
plotted across \(n_{\mathrm{obs}}\). Panel/label/reference-band
conventions match Figures 1-2, except the grouping is by variant rather
than by SYDNI-vs-comparator: each panel's gray band and dashed line trace
the min-max range and mean of the *other variant* at each
\(n_{\mathrm{obs}}\), so the flat-vs-ma gap that is this study's central
finding is visible directly within the figure rather than requiring a
cross-reference to a separate panel. The `flat` panel sitting well above
its own reference band (the `ma` range) on every metric, and the `ma`
panel sitting at or below its reference band (the `flat` range), is the
same result reported numerically above.

**Scenario 3: indicator-level data exists and item-level output is
needed.** This is where the calibration mechanism's behavior can be judged
on its own terms, ranking the five `ma` engines by mean rank across Type I
error, false discovery rate, power recovery, within-scale correlation
error, between-construct sign-reversal rate, and runtime:

\begin{center}
\begin{tabular}{lrrrrrrr}
\hline
Engine & Mean rank & Type I & FDR & Power & Within-scale error & Between sign-rev.\ & sec/condition \\
\hline
\texttt{sydni\_factor\_ma}    & 1.58 & 0.028 & 0.015 & 0.022 & 0.031 & 0.041 & 1.01 \\
\texttt{sydni\_graphical\_ma} & 2.08 & 0.033 & 0.018 & 0.039 & 0.034 & 0.041 & 0.42 \\
\texttt{sydni\_copula\_ma}    & 3.17 & 0.039 & 0.021 & 0.028 & 0.050 & 0.042 & 0.46 \\
\texttt{sydni\_bart\_ma}      & 3.67 & 0.039 & 0.021 & 0.061 & 0.117 & 0.027 & 1.39 \\
\texttt{sydni\_brm\_ma}       & 4.50 & 0.050 & 0.028 & 0.050 & 0.042 & 0.044 & 585.9 \\
\hline
\end{tabular}
\end{center}

Two findings here revise this family's Section 8.6-derived positioning,
scoped strictly to Scenario 3. First, **`sydni_brm`'s measurement-aware
path is not the all-arounder its manifest-variable behavior would
suggest.** Where Section 8.6 found `sydni_brm` never worst on any axis,
here it finishes last by mean rank and is worst or tied-worst on four of
six axes - Type I error, false discovery rate, between-construct
sign-reversal, and runtime, the last at 400-1,400x the compute of the
other four engines - with only power recovery (second of five) as a clear
strength. This does not touch Section 8.6's manifest-variable conclusion
about `sydni_brm`; it is confined to the measurement-aware use case.
Second, **`sydni_bart`'s niche narrows for this scenario without being
eliminated.** It still posts the best power recovery (0.061) and the best
between-construct sign-reversal rate (0.027) of the five, but its
within-scale correlation error (0.117) is 3-4x every other `ma` engine's
(0.031-0.050) - the one engine whose measurement-aware path does not
clearly protect the within-construct reliability the mechanism is
specifically designed to protect. Where Section 8.6 positions
`sydni_bart` as a general near-free `sydni_brm` stand-in, for
measurement-aware use it is better read as a specialist for when power or
between-construct fidelity matters more than item-level reliability, not
as a default. `sydni_factor_ma` and `sydni_graphical_ma` emerge as the
recommended defaults for this scenario instead - both rank in the top
half of the family on every single axis, at 0.42-1.01 seconds per
condition; `sydni_factor_ma` trades the most power for the strongest
safety numbers, `sydni_graphical_ma` is the faster, more balanced
alternative of the two.

Combining all three scenarios gives a scenario-dependent recommendation
rather than one global engine ranking:

\begin{center}
\begin{tabular}{p{4.2cm}p{5.3cm}p{2.8cm}}
\hline
User need & Recommended engine(s) & Evidence source \\
\hline
Manifest/composite predictors, no indicator data & \texttt{sydni\_brm} (reference); \texttt{sydni\_copula} (fast primary); \texttt{sydni\_bart}/\texttt{sydni\_graphical} (cheap secondary); \texttt{sydni\_factor} (safety extreme) & Section 8.6 \\
Indicator-level data, composite-score output only & Supply \texttt{measurement\_blocks} regardless of engine & Section 8.7 (flat vs.\ ma) \\
Indicator-level data, item-level/SEM output needed & \texttt{sydni\_factor\_ma}/\texttt{sydni\_graphical\_ma} (defaults); \texttt{sydni\_bart\_ma} (power/between-construct specialist); \texttt{sydni\_brm\_ma} (viable but expensive, no demonstrated payoff) & Section 8.7 (ma arm) \\
\hline
\end{tabular}
\end{center}

**Sample-size dependence** mirrors Section 8.6's pattern on this different
DGP: power recovery for every `ma` engine shrinks toward zero as
\(n_{\mathrm{obs}}\) grows (e.g. `sydni_brm_ma`: 0.117 to 0.017 from
\(n_{\mathrm{obs}} = 75\) to 300), and each `ma` engine's Type I error
converges toward the observed-only baseline's own rate by
\(n_{\mathrm{obs}} = 300\), consistent with SYDNI's positioning as a
small/underpowered-sample tool rather than a general-purpose one.

**Open questions carried forward from this pilot.** Why `sydni_bart`'s
measurement-aware path does not protect within-scale reliability the way
the other four engines' do is not diagnosed here - a documented open
question alongside Section 8.6's still-open BART nonlinearity question,
not a resolved one. `sydni_brm_ma`'s runtime cost was traced during this
pilot's setup to a `brms`/Stan mechanics issue rather than a modeling
choice: a per-replicate, data-estimated intercept-prior value is embedded
as a literal directly in the generated Stan program text, forcing a
near-total compilation-cache miss on almost every call regardless of
predictor selection; this is noted as a known cost driver rather than
something correctable from the package's benchmarking layer. Finally,
unlike Section 8.6's seven-scenario design, every finding in this section
comes from one latent-structure DGP, and cannot speak to how the flat-vs-ma
gap or the ma-engine ranking would generalize to a different measurement
model.

## 9. Empirical illustration

The package's included `grit` dataset (500 respondents' Big Five
personality items and Duckworth et al.'s (2007) Grit Scale, from
OpenPsychometrics) provides a natural empirical example, because it lets
Sections 8.6-8.7's simulation-based findings be checked against a real
psychometric instrument rather than only a designed data-generating
process. This section reports a single-run illustration (Study 3), not a
replicated benchmark: real data has no known ground truth, so the metrics
below are structure-retention and downstream-model-stability comparisons
rather than the Type I error/power-recovery metrics Sections 8.6-8.7 use
against known truth. Full results, reproducing scripts, and additional
detail beyond this section's summary are reported in
`results/case_studies/grit_illustration/FINDINGS.md` in the package
repository.

**Design.** 447 complete cases (of 500) were used, with Conscientiousness
excluded from the predictor set: its correlation with the Grit composite
(r = 0.61) is more than double any other Big Five trait's (Neuroticism r
= -0.40, Agreeableness r = 0.24, Extraversion r = 0.19, Openness r =
0.13), consistent with a documented construct-redundancy critique in the
grit literature (Credé, Tynan, & Harms, 2017) rather than a
statistical concern - Conscientiousness's correlation with the remaining
four traits is modest (|r| <= 0.27). The four remaining Big Five
batteries (10 items each) served as predictor item blocks. Grit was
modeled using the hierarchical `measurement_blocks` composite-list
feature from Section 4.1 for the first time in this paper: Consistency of
Interest (6 items) and Perseverance of Effort (6 items) as lower-order
leaves, combined into one higher-order Grit structural family. This is the
`list(lower = ..., higher = ...)` composite-list form of
`measurement_blocks` mentioned in Section 4.1: the four Big Five batteries
remain flat leaf/final families, while Consistency of Interest and
Perseverance of Effort are declared as leaves under a single higher-order
`Grit` family:

```r
measurement_blocks <- list(
  lower = list(
    E = paste0("e", 1:10), N = paste0("n", 1:10),
    A = paste0("a", 1:10), O = paste0("o", 1:10),
    CI = c("gs2", "gs3", "gs5", "gs7", "gs8", "gs11"),
    PE = c("gs1", "gs4", "gs6", "gs9", "gs10", "gs12")
  ),
  higher = list(Grit = c("CI", "PE"))
)

ma_out <- sydni(
  data = grit_complete,
  variables = c(paste0("e", 1:10), paste0("n", 1:10), paste0("a", 1:10),
                paste0("o", 1:10), paste0("gs", 1:12)),
  engine = "copula",
  measurement_blocks = measurement_blocks,
  between_block_shrink = 0.5,
  seed = 123
)
```

(The actual illustration also included a `class` covariate - a median
split on the Grit composite, added only to satisfy `SMOTE`/`ADASYN`/
`synthpop`'s hardcoded binary-target dependency for the off-label
comparators; it plays no role in the `measurement_blocks` mechanism shown
above and is omitted here for clarity.)

All 17 methods from Section 8.7 (the 7 off-label comparators plus each
SYDNI engine as `_flat` and `_ma`) were run once each this way (`n_syn`
equal to `n_obs`), with 0 failures. The downstream model,
\(\mathrm{Grit} \sim E + N + A + O\), found Neuroticism (b = -0.30,
p < .001) and Agreeableness (b = 0.16, p < .001) significant in the
observed data, with Extraversion (b = 0.005, p = .88) and Openness
(b = 0.076, p = .10) not significant.

**Downstream-model stability.** Every method that flagged a term as
newly significant after augmentation flagged it on Openness - the one
borderline predictor in the observed data (p = .10), not a term that was
strongly non-significant to begin with. The *size* of the change split
cleanly by mode: `sydni_bart_ma` and `sydni_factor_ma` moved Openness's
point estimate almost imperceptibly (0.076 to 0.071 and 0.076 to 0.076)
while its p-value dropped from a smaller standard error - the pattern one
would expect from recovering power on an underpowered real effect, though
real data with no known ground truth cannot distinguish that from a
smaller, coincidentally power-favorable draw; Sections 8.6-8.7's
simulations, which do have known truth, are what actually establish that
distinction. `_flat`-mode engines moved the estimate further
(`sydni_brm_flat`: 0.076 to 0.109, +43%). **`sydni_graphical_flat`
produced the one result worth real caution**: it turned an Extraversion
effect that was essentially flat in the observed data (0.005, p = .88)
into a borderline "finding" (0.058, p = .048) - `sydni_graphical_ma`
shows nothing of the kind, a concrete real-data demonstration of Section
8.7's "always supply `measurement_blocks` when indicator-level data
exists" recommendation.

**Structure retention replicates Section 8.7's core finding on real
data.** For `copula`, `graphical`, `factor`, and `brm`, `_ma` mode cuts
within-scale distortion sharply relative to `_flat` (e.g. `graphical`:
within-scale correlation error 0.138 to 0.032) and flips between-construct
correlation change from *inflating* under `_flat` (three of the four
engines inflated between-construct correlations on this real dataset -
not merely a simulated risk) to properly *shrinking* under `_ma`.
**`sydni_bart` is again the exception, more starkly than in Section
8.7:** both within-scale metrics get *worse* going from `_flat` to `_ma`
(within-scale correlation error 0.080 to 0.097) - the only engine where
`_ma` moves the wrong direction on real data, and unlike Section 8.7
(where `sydni_bart_ma` was merely the weakest of five `_ma` engines),
here `_ma` is worse than the same engine's own `_flat` version outright.
The simulated and single-run real-data results are consistent in showing
weaker within-scale retention for `sydni_bart_ma`, but they do not establish
a general engine property: the simulation uses one measurement model and the
empirical illustration has no replication or known population truth.

**Comparison against the off-label incumbents.** `ADASYN`, `bootstrap`,
and `SMOTE` post lower raw within-scale correlation error than any `_ma`
engine - the same row-copying/light-interpolation artifact documented in
Sections 8.6-8.7, not genuine structure modeling. `synthpop_sample`
performs worst of all 17 methods on within-scale retention (correlation
error 0.187, worse than every SYDNI variant including `_flat`), since its
marginal-only synthesis has no mechanism to preserve cross-item
dependence at all. None of the off-label comparators moved a downstream
term that was strongly non-significant in the observed data to
significance in this run, but this single,
unreplicated illustration cannot itself establish a false-discovery rate
the way Sections 8.6-8.7's replicated designs do; it is consistent with,
but not independent confirmation of, those sections' finding that methods
lacking any shrinkage mechanism carry a false-positive cost that
requires many replications under known truth to detect reliably.

**Runtime** matched the established cost profile: `sydni_brm_flat` (52
separate per-item model fits) took approximately 73 minutes on this
illustration; `sydni_brm_ma` (5 family-level models) took approximately 5
minutes - the largest `_flat`-vs-`_ma` runtime gap of any engine, and a
further real-data demonstration of why `brm` remains the most expensive
of the five engines regardless of mode.

This illustration should be read as a worked example that happens to
corroborate and sharpen Sections 8.6-8.7's simulation-based findings, not
as an independently powered confirmatory result - a single real dataset
with no replication cannot carry the same evidential weight as either
simulation study's grid. Disagreement across engines here is informative
in the same way Sections 8.6-8.7 treat it: it reveals which conclusions
depend on a particular engine's generative assumptions rather than
holding regardless of synthesis mechanism.

## 10. Reproducibility and software availability

The draft corresponds to `sydni` package version `0.2.0`. The package depends
on R version `>= 4.4.0`. Core imports include `dplyr`, `purrr`, `truncnorm`,
`copula`, `corpcor`, `infotheo`, and `progress`.

Optional engine dependencies include:

- `brms` for `engine = "brm"`;
- `dbarts` for `engine = "bart"`;
- `psych` for `engine = "factor"`;
- `huge` for `engine = "graphical"`;
- `cmdstanr` for BRM models using the `cmdstanr` backend.

`brms` is a `Suggests` dependency rather than a hard import, so the package
can be installed and used with the other four engines without installing
`brms` or a Stan toolchain.

The package includes a benchmark scaffold under `inst/benchmarks/` and
methodology vignettes under `vignettes/`. A submission version of this paper
should include reproducible benchmark scripts, session information, and all
simulation settings needed to regenerate tables and figures.

## 11. Limitations

SYDNI has several important limitations:

1. It cannot create empirical information absent from the observed sample.
2. It does not provide formal privacy or disclosure protection.
3. Very small subgroup cells can lead to unstable subgroup structure.
4. Ordinal variables require careful validation, especially when treated as
   numeric approximations by conditional engines.
5. Copula and graphical engines may smooth or miss tail behavior.
6. Factor synthesis is domain-specific and best suited to item batteries or
   related measurement settings.
7. BRM and BART engines depend on the adequacy of their target-wise conditional
   model assumptions. Pilot benchmarking (Section 8.6), including a direct
   follow-up test using an oracle analysis model that explicitly credits
   correct nonlinear/interaction recovery, did not find BART's theoretical
   nonlinear/interaction flexibility translating into a measurable
   inferential-safety advantage over the other engines on manifest
   variables; this is treated as a demonstrated non-finding rather than a
   settled limitation, since it does not rule out an advantage on data
   this benchmark's scenarios do not stress-test. In measurement-aware mode
   specifically, a second pilot (Section 8.7) found BART's within-scale
   item-reliability retention 3-4x worse than the other four engines',
   even though it led the family on power recovery and between-construct
   fidelity there - a real, engine-specific tradeoff rather than a
   universal weakness. A real-data illustration (Section 9) reproduced
   this weakness independently, and more starkly: within-scale retention
    got worse, not merely weaker-than-average, when moving from BART's
    flat to measurement-aware mode. The simulated and single-run real-data
    results are consistent with this weakness, but one measurement-model DGP
    and one unreplicated empirical illustration do not establish a general
    engine property.
8. The computational cost differs substantially across engines, especially for
   BRM models requiring Stan backends; in measurement-aware mode this cost
   is further amplified (Section 8.7) by `brms` embedding a per-replicate,
   data-estimated prior value as a literal in the generated Stan program,
   which forces a near-total compilation-cache miss on almost every call
   and was not found to correspond to a measurable quality advantage over
   the other engines' measurement-aware paths in that pilot.
9. In BRM's measurement-aware mode, structural-family-score calibration is
   applied before item-level donor reconstruction, so final item-level
   composite correlations can differ slightly from the calibration target,
   unlike the other four engines.
10. Retention-aware group balancing currently supports only
    `balance_target = "largest_group"`.

These limitations are consistent with the intended use of SYDNI as a
robustness workflow, not a universal synthetic-data generator.

## 12. Discussion

SYDNI frames synthetic augmentation as a conservative sensitivity-analysis
tool. The key question is not whether synthetic rows should be treated as
additional observations, but whether analytic conclusions are stable when the
observed sample is extended through controlled perturbation and
structure-preserving synthesis.

The five-engine design allows analysts to ask different robustness questions.
If a conclusion is stable under BRM and BART, it may be less sensitive to
linear versus nonlinear conditional modeling choices - though Section 8.6's
direct test of this premise found no measurable empirical difference between
the two on the scenario constructed to separate them, so this framing should
currently be read as the design intent rather than a demonstrated
distinction. If it is stable under
copula and graphical engines, it may be less sensitive to dense versus sparse
dependence representations. If a psychometric conclusion is stable under the
factor engine, it may be less sensitive to latent-structure perturbation.
Conversely, disagreement across engines can be substantively useful: it may
show that an apparent result depends on a particular view of the data-generating
structure.

Two pilot-scale simulation benchmarks have now been completed. The
manifest-variable pilot (Section 8.6) supports the
conservative-inferential-safety premise above relative to off-label
resampling comparators. It also directly tested, rather than left open,
whether BART shows a nonlinear/interaction advantage: an oracle-form
analysis model built specifically to credit correct recovery of the
`nonlinear` scenario's true threshold and interaction terms still found no
measurable BART advantage over the other four engines. This is treated as
a demonstrated non-finding under this benchmark's specific DGP and
analysis target, not as a general claim that BART's theoretical
nonlinear/interaction flexibility never matters - it narrows what future
testing on different nonlinear data shapes would need to show, rather than
closing the question outright. The measurement-aware calibration pilot
(Section 8.7) directly tests the mechanism described in Section 4.1 for
the first time and shows it working as designed for four of the five
engines, while revising the engine-recommendation picture in a
scenario-dependent way: `sydni_brm`'s manifest-variable positioning as the
balanced reference engine does not carry over to measurement-aware use,
where it is dominated on cost and on most safety/structure axes by
`sydni_factor` and `sydni_graphical`, and `sydni_bart`'s role narrows from
a general `sydni_brm` stand-in to a specialist for power- and
between-construct-priority cases. Neither result should be read as an
overall winner across both pilots - which engine to recommend depends on
whether the analyst has manifest/composite variables only, indicator-level
data used only for composite scores, or indicator-level data analyzed
directly, as summarized in Section 8.7.

A single-run empirical illustration on the package's real `grit` dataset
(Section 9) checked both pilots' central findings against a genuine
psychometric instrument rather than only designed data-generating
processes, and corroborated them: measurement-aware calibration protected
within-scale structure and shrank between-construct associations for
four of the five engines exactly as Section 8.7 found, one off-label
comparator run without calibration turned a near-null effect into a
borderline "finding" in a way its measurement-aware SYDNI counterpart did
not, and BART's within-scale weakness reproduced independently - more
starkly than in simulation. This single real dataset cannot carry the
same evidential weight as either pilot's replicated grid, but the
alignment between simulated and real-data results in the same direction
is itself informative: none of Sections 8.6-8.7's central conclusions
depended on an artifact specific to their designed scenarios.

Future work should prioritize testing generalizability over adding raw
replications: both pilots' large directional findings (SYDNI versus
naive resampling, the measurement-aware mechanism's effect,
BART's within-scale weakness) are already well-supported at 10
replications per condition, since they reflect Monte Carlo error that
more reps would shrink, not the dominant source of uncertainty here.
The dominant source is external validity - how well any single simulated
design matches the range of real data analysts actually encounter - which
more replications of the *same* design cannot address. The more valuable
next step is therefore testing BART's nonlinear/interaction premise and
the measurement-aware mechanism against additional, structurally
different data-generating designs and additional real datasets beyond
Section 9's single illustration, not scaling the existing designs to
more replications. Remaining priorities include adding further published
examples, improving diagnostics for synthetic-data quality, exploring
missing-data integration, and adding optional disclosure-risk diagnostics
for users who need privacy-oriented assessments.

## Appendix A: Engine-specific pseudocode

Section 5 describes each engine's design and role at a level suited to
reading the main text linearly. This appendix gives full step-by-step
pseudocode for all five engines, for readers who want the exact generation
sequence.

### A.1 BRM engine

```text
Algorithm 2: BRM engine

Input: observed data D, variables V, covariates C, settings theta

1. Draw bootstrap sample D_boot of size n_aug.
2. Add noise to synthesized columns in D_boot -> D_feed.
3. For each target variable v in V:
   a. Select predictors Z_v from V \ v.
   b. Append covariates C to the predictor set.
   c. Choose or infer outcome family F_v.
   d. Build augmentation-linked shrinkage prior for slopes.
   e. Build intercept prior from the observed target distribution.
   f. Fit brms model M_v on D_feed.
   g. Generate posterior predictive draws for D_feed.
   h. Reconstruct output to original variable type and support.
4. Combine synthetic target columns.
5. Append covariates from D_feed.
6. Return models and synthetic data.
```

### A.2 BART engine

```text
Algorithm 3: BART engine

Input: observed data D, variables V, covariates C, settings theta

1. Draw bootstrap sample D_boot of size n_aug.
2. Add noise to synthesized columns in D_boot -> D_feed.
3. Classify targets as continuous, binary, or ordinal-like.
4. For each target variable v in V:
   a. Select predictors Z_v from V \ v.
   b. Build a numeric design matrix from Z_v and C.
   c. Compute augmentation-linked BART regularization.
   d. Fit BART model M_v on D_feed.
   e. Draw posterior predictions for feeder rows.
   f. Threshold binary targets or reconstruct ordinal-like targets.
5. Combine synthetic target columns.
6. Append covariates from D_feed.
7. Return models and synthetic data.
```

### A.3 Copula engine

```text
Algorithm 4: Copula engine

Input: observed data D, variables V, covariates C, settings theta

1. Partition D into subgroups by covariate combinations.
2. Compute global shrinkage correlation matrix R_global.
3. For each subgroup g:
   a. Determine target feeder size.
   b. Bootstrap additional subgroup rows if needed.
   c. Add noise to synthesized variables.
   d. Pool noisy rows with observed subgroup rows.
   e. Compute subgroup shrinkage correlation matrix R_g.
   f. Blend R_g with R_global and stabilize the result.
   g. Simulate Gaussian copula uniforms U_g.
   h. Map U_g through empirical subgroup margins.
   i. Append subgroup covariates.
4. Combine subgroup outputs.
5. Return subgroup models and synthetic data.
```

### A.4 Factor engine

```text
Algorithm 5: Factor engine

Input: observed data D, variables V, covariates C, settings theta

1. Partition D into subgroups by covariate combinations.
2. Build noisy feeder data within each subgroup.
3. For each subgroup g:
   a. Identify active non-constant subgroup items A_g.
   b. Estimate a subgroup correlation matrix.
   c. Fit a factor model on A_g when viable.
   d. Shrink loadings, factor correlations, and uniquenesses.
   e. Generate latent factor draws and residuals.
   f. Construct latent item values.
   g. Map active items back through empirical subgroup margins.
   h. Reconstruct inactive items from observed subgroup support.
   i. If factor fitting fails, fall back to empirical subgroup reconstruction.
4. Combine subgroup outputs.
5. Return subgroup models and synthetic data.
```

### A.5 Graphical engine

```text
Algorithm 6: Graphical engine

Input: observed data D, variables V, covariates C, settings theta

1. Partition D into subgroups by covariate combinations.
2. Build noisy feeder data within each subgroup.
3. For each subgroup g:
   a. Transform subgroup feeder data to latent Gaussian space.
   b. Identify active variables with non-zero latent variance.
   c. Estimate sparse precision structure on active variables.
   d. Convert precision matrix to stabilized correlation matrix.
   e. If sparse estimation fails, fall back to stabilized latent correlation.
   f. Simulate Gaussian latent draws.
   g. Transform latent draws to uniforms.
   h. Map uniforms back through empirical subgroup margins.
   i. Append subgroup covariates.
4. Combine subgroup outputs.
5. Return subgroup models and synthetic data.
```

## 13. References

References below are organized by topic. Entries for R packages give the
CRAN landing page URL (stable across releases); the exact version number
and citation year should be regenerated via `citation("<package>")`
immediately before submission, since both drift with each CRAN release
and are not pinned here.

### Methodological and psychometric literature

Chawla, N. V., Bowyer, K. W., Hall, L. O., & Kegelmeyer, W. P. (2002).
SMOTE: Synthetic minority over-sampling technique. *Journal of Artificial
Intelligence Research*, 16, 321-357.

Chipman, H. A., George, E. I., & McCulloch, R. E. (2010). BART: Bayesian
additive regression trees. *The Annals of Applied Statistics*, 4(1),
266-298.

Credé, M., Tynan, M. C., & Harms, P. D. (2017). Much ado about grit: A
meta-analytic synthesis of the grit literature. *Journal of Personality
and Social Psychology*, 113(3), 492-511.

Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
tests. *Psychometrika*, 16(3), 297-334.

Duckworth, A. L., Peterson, C., Matthews, M. D., & Kelly, D. R. (2007).
Grit: Perseverance and passion for long-term goals. *Journal of
Personality and Social Psychology*, 92(6), 1087-1101.

Efron, B. (1979). Bootstrap methods: Another look at the jackknife. *The
Annals of Statistics*, 7(1), 1-26.

Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
covariance estimation with the graphical lasso. *Biostatistics*, 9(3),
432-441.

He, H., Bai, Y., Garcia, E. A., & Li, S. (2008). ADASYN: Adaptive
synthetic sampling approach for imbalanced learning. In *2008 IEEE
International Joint Conference on Neural Networks (IEEE World Congress on
Computational Intelligence)* (pp. 1322-1328). IEEE.

Iman, R. L., & Conover, W. J. (1982). A distribution-free approach to
inducing rank correlation among input variables. *Communications in
Statistics - Simulation and Computation*, 11(3), 311-334.

John, O. P., & Srivastava, S. (1999). The Big Five trait taxonomy:
History, measurement, and theoretical perspectives. In L. A. Pervin &
O. P. John (Eds.), *Handbook of Personality: Theory and Research* (2nd
ed., pp. 102-138). Guilford Press.

Nelsen, R. B. (2006). *An Introduction to Copulas* (2nd ed.). Springer.

Olsson, U. (1979). Maximum likelihood estimation of the polychoric
correlation coefficient. *Psychometrika*, 44(4), 443-460.

Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys*.
Wiley.

van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate
imputation by chained equations in R. *Journal of Statistical Software*,
45(3), 1-67.

Zhao, T., Liu, H., Roeder, K., Lafferty, J., & Wasserman, L. (2012). The
`huge` package for high-dimensional undirected graph estimation in R.
*Journal of Machine Learning Research*, 13, 1059-1062.

### Software

Bürkner, P.-C. (2017). brms: An R package for Bayesian multilevel models
using Stan. *Journal of Statistical Software*, 80(1), 1-28.

Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B.,
Betancourt, M., Brubaker, M., Guo, J., Li, P., & Riddell, A. (2017).
Stan: A probabilistic programming language. *Journal of Statistical
Software*, 76(1), 1-32.

Csárdi, G., & FitzJohn, R. progress: Terminal Progress Bars. R package.
https://CRAN.R-project.org/package=progress

Dorie, V. dbarts: Discrete Bayesian Additive Regression Trees Sampler. R
package. https://CRAN.R-project.org/package=dbarts

Gabry, J., Češnovar, R., & Johnson, A. cmdstanr: R Interface to
'CmdStan'. https://mc-stan.org/cmdstanr/

Henry, L., & Wickham, H. purrr: Functional Programming Tools. R package.
https://CRAN.R-project.org/package=purrr

Hofert, M., Kojadinovic, I., Maechler, M., & Yan, J. copula: Multivariate
Dependence with Copulas. R package.
https://CRAN.R-project.org/package=copula

Mersmann, O., Trautmann, H., Steuer, D., & Bornkamp, B. truncnorm:
Truncated Normal Distribution. R package.
https://CRAN.R-project.org/package=truncnorm

Meyer, P. E. infotheo: Information-Theoretic Measures. R package.
https://CRAN.R-project.org/package=infotheo

Nowok, B., Raab, G. M., & Dibben, C. (2016). synthpop: Bespoke creation
of synthetic data in R. *Journal of Statistical Software*, 74(11), 1-26.

R Core Team. *R: A Language and Environment for Statistical Computing*.
R Foundation for Statistical Computing, Vienna, Austria.
https://www.R-project.org/

Revelle, W. psych: Procedures for Psychological, Psychometric, and
Personality Research. Northwestern University, Evanston, Illinois. R
package. https://CRAN.R-project.org/package=psych

Schafer, J., Opgen-Rhein, R., Zuber, V., Ahdesmaki, M., Silva, A. P. D.,
& Strimmer, K. corpcor: Efficient Estimation of Covariance and (Partial)
Correlation. R package. https://CRAN.R-project.org/package=corpcor

Siriseriwan, W. smotefamily: A Collection of Oversampling Techniques for
Class Imbalance Problem Based on SMOTE. R package.
https://CRAN.R-project.org/package=smotefamily

Stan Development Team. RStan: The R Interface to Stan.
https://mc-stan.org/rstan/

Wickham, H., François, R., Henry, L., & Müller, K. dplyr: A Grammar of
Data Manipulation. R package. https://CRAN.R-project.org/package=dplyr

### Data

The grit case study (Section 9) uses the package's bundled `grit`
dataset, itself drawn from Big Five personality items and Duckworth et
al.'s (2007) Grit Scale; see Duckworth et al. (2007) and John and
Srivastava (1999) above for the underlying instruments.
