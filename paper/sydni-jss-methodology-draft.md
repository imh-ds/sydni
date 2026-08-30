# sydni: Synthetic Data Augmentation for Sensitivity Analysis in R

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
(SYDNI), an R package for sensitivity analysis through controlled synthetic
perturbation. It constructs noisy bootstrap feeder data, applies one of five
swappable synthesis engines, and uses augmentation-linked shrinkage so that
synthetic rows are not treated as ordinary additional empirical evidence.
The intended use is a pre-specified comparison of the original analysis with
analyses under plausible synthesis choices: agreement can support a claim of
stability, whereas disagreement identifies assumptions to investigate. The
package supports conditional, dependence-first, and latent-structure engines,
with optional measurement-aware calibration and retention-aware augmentation.
This article describes the workflow, software interface, validation and
reporting protocol, and pilot validation studies. A single psychometric-data
example illustrates the workflow; it is not treated as confirmatory evidence.

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
structure-preserving synthesis engine. Here, “conservative” is a design goal:
the workflow is intended to make augmentation-induced changes visible and
auditable, not to establish a formal guarantee about Type I error or effect
recovery in every application.

This use differs from multiple imputation, whose aim is inference with missing
data (Rubin, 1987; van Buuren & Groothuis-Oudshoorn, 2011), and from synthetic
microdata systems such as `synthpop`, whose central aim is data release
(Nowok, Raab, & Dibben, 2016). It also differs from class-imbalance procedures
such as SMOTE and ADASYN, which target predictive classification rather than
the stability of an observational estimand (Chawla et al., 2002; He et al.,
2008). These methods can be informative boundary comparisons, but they are not
interchangeable inferential tools.

The following table positions SYDNI against the closest adjacent tools. We are not aware
of an open-source R package whose primary purpose is synthetic augmentation as
a non-confirmatory sensitivity workflow. The absence of a direct comparator
does not make adjacent methods irrelevant: they are useful baselines for
showing how tools designed for a different task behave when analysts use them
for augmentation.

| Tool family | Primary intended use | Relation to SYDNI | Role in this paper |
|---|---|---|---|
| SYDNI | Sensitivity analysis under controlled synthetic perturbation | The focal workflow; original-data inference remains primary | Evaluated across engines and design settings |
| Multiple imputation | Inference with missing data | Addresses missingness rather than augmentation sensitivity | Positioned conceptually; not a direct benchmark |
| Synthetic microdata (`synthpop`) | Privacy-oriented or disclosure-controlled data release | Can generate synthetic records, but has a different inferential target | Adjacent baseline, evaluated outside its intended use |
| Bootstrap resampling | Sampling-variation assessment | Resamples observed data without SYDNI's synthesis and shrinkage mechanisms | Minimal augmentation baseline |
| SMOTE and ADASYN | Class-imbalance learning | Upsample minority classes for predictive tasks, not observational estimands | Off-label stress-test baselines |
| Proprietary augmentation services | Product-specific data generation or augmentation | Algorithms and reproducibility conditions may be unavailable for scholarly comparison | Not benchmarked |

The comparisons in this paper therefore assess adjacent or off-label baselines,
not head-to-head substitutes for the same inferential task. They should not be
read as claims that SYDNI is superior for privacy-preserving release,
classification, or missing-data inference.

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
increases evidential sample size. Its contribution is a perturbation framework
for asking whether analytic conclusions are fragile to
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
single-run empirical illustration on a real psychometric dataset (Section 9).

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
3. Construct a feeder dataset by bootstrap resampling observed rows (Efron,
   1979).
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
noisy subgroup feeder data, estimates subgroup Gaussian copula correlations
(Nelsen, 2006),
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
| `brm` | Conditional regression | Native predictors | Explicit parametric models and mixed-outcome support | Computational cost and model specification |
| `bart` | Nonlinear conditional regression | Native predictors | Flexible conditional functions and interactions | Ordinal targets use an approximation path |
| `copula` | Dense joint dependence | Subgrouping | Empirical margins with a joint dependence model | Tail behavior may be smoothed |
| `factor` | Latent measurement structure | Subgrouping | Psychometric item batteries | Sparse subgroups can be fragile |
| `graphical` | Sparse direct dependence | Subgrouping | Network-like conditional structure | Sparsity must be scientifically meaningful |

All five engines additionally accept the measurement-aware calibration and
retention-aware augmentation controls described in Sections 4.1 and 4.2. The
underlying algorithmic target in the table above (conditional regression,
dense dependence, and so on) does not change; measurement-aware mode adds a
between-construct shrinkage and calibration layer on top of each engine's
existing generative mechanism. Comparative empirical results are reported in
Section 8; they are pilot-specific and should not be read as a global engine
ranking. Engine choice should follow the data representation, scientific
dependence assumptions, and computational budget.
Full step-by-step pseudocode for all five engines is given in Appendix A.

## 6. Software interface

The public package interface centers on `sydni()`. A typical call requires a
dataset, the variables to synthesize, an engine, an augmentation size, and a
seed:

```r
sydni(
  data,
  variables,
  covariates = NULL,
  engine = "copula",
  n_aug,
  noise_method = "truncnormal",
  measurement_blocks = NULL,
  between_block_shrink = NULL,
  seed
)
```

The complete argument reference, including noise, predictor, retention,
balancing, and engine-specific controls, is provided in `?sydni` and the
package vignette. The returned list contains `data` (the generated rows),
`models` (engine-specific fitted objects), and `meta` (generation settings and
diagnostics). A manuscript workflow should report the settings in `meta`, not
only a shortened example call.

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

str(result, max.level = 1)
```

```
List of 3
 $ models:List of 1
 $ data  :'data.frame': 100 obs. of 22 variables
 $ meta  :List of generation settings and diagnostics
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

The following table gives a reporting and interpretation protocol. Before generating data,
analysts should predeclare the engines, settings, seeds, and the substantively
meaningful tolerance used to judge effect-size changes. No universal numerical
tolerance is implied: it depends on the estimand and scientific context.

| Original analysis | Sensitivity results across the predeclared set | Interpretation and reporting |
|---|---|---|
| Result is present | Direction, practical magnitude, interval interpretation, and decision are similar | Report stability under the examined perturbations; do not describe this as stronger evidence |
| Result is absent | Result appears only after one or more augmentations | Report an exploratory sensitivity signal, not a confirmed finding |
| Result is present | Direction, practical magnitude, interval interpretation, or decision changes materially | Report fragility and investigate the implicated data or modeling assumption |
| Any original result | Results differ materially across engines, settings, or seeds | Report model dependence; do not select the most favorable augmented run |

Recommended reporting language:

> SYDNI analyses were treated as secondary robustness checks alongside the
> original-sample analyses. Synthetic augmentation was not interpreted as
> additional empirical evidence. We predeclared the sensitivity set and effect
> tolerance, reported the original analysis first, and treated augmentation-only
> effects or disagreement across runs as exploratory sensitivity signals.

Manuscripts using SYDNI should report the full generation settings, including
`engine`, `n_aug` or `prop_aug`, `noise_method`, noise controls, predictor
selection settings, shrinkage settings, backend settings, and seed. They should
also report the original analysis, the complete predeclared sensitivity set,
effect and interval comparisons, and the diagnostics used to interpret any
disagreement.

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
truth. Both benchmarks (Sections 8.6-8.7) used:

\[
n_{\mathrm{obs}} \in \{75, 150, 300\},
\qquad
n_{\mathrm{syn}} / n_{\mathrm{obs}} \in \{0.5, 1, 2\}.
\]

Both began as 10-replication pilots, sufficient for a compact comparison
but not for precise Monte Carlo estimates or a stable fine-grained ranking
among close methods; point summaries from a 10-replication run should be
read as directional evidence within these data-generating designs rather
than resolved differences. Both have since been confirmed at much higher
replication (100 per condition for every method except `sydni_brm`, 50
for `sydni_brm` on a narrower focused grid, in each study - see Sections
8.6 and 8.7 for the precision this achieves).

Section 8.6's manifest-variable benchmark covers seven simulation scenarios,
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
- SMOTE and ADASYN from `smotefamily`, used as off-label class-imbalance
  baselines with a binary balancing target;
- `synthpop`, in both a sequential/CART configuration and a conservative
  sample-based configuration, as an adjacent synthetic-microdata baseline;
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
- Polychoric correlation error for Likert items (Olsson, 1979).
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

Section 8.6 reports confirmation-scale numbers for Tables 1-3 and Figures
1-3 on manifest variables; Section 8.7 reports a companion confirmation
study isolating the measurement-aware calibration mechanism on a
latent-variable item-battery design, including Figure 4; Section 9 reports
a single-run empirical illustration. Section 12 discusses the need for
additional data-generating designs beyond the two exercised here to assess
generalizability.

### 8.6 Confirmed results

An initial 10-replication pilot of the design above (all 7 scenarios from
Section 8.1, the full \(n_{\mathrm{obs}} \in \{75, 150, 300\}\) and
augmentation-ratio grid, 12 methods spanning every category in Section 8.2
except a tabular GAN comparator, 7,560 total conditions) surfaced the
pattern reported below, but 10 replications cannot support precise Monte
Carlo estimates or a stable fine-grained ranking among close methods
(Section 8.1). A confirmation run was therefore completed at substantially
higher replication and is reported as the primary result throughout this
section: 100 replications per condition for every method except
`sydni_brm` (69,300 conditions, the same full grid as the pilot), and 50
replications on a predeclared focused subset for `sydni_brm` alone - 5 of
the 7 scenarios (`null_only`, `sparse_effects`, `correlated_predictors`,
`moderation_imbalance`, `small_cell`), \(n_{\mathrm{obs}} \in \{75, 150\}\),
\(n_{\mathrm{syn}}/n_{\mathrm{obs}} \in \{1, 2\}\), 1,000 conditions.
`sydni_brm` was confirmed on a narrower grid because it costs roughly
three orders of magnitude more compute per condition than every other
method (Figure 3); at 50 replications its design-stratified 95% Monte
Carlo interval (MCI) half-width - computed per design cell and combined
across the predeclared grid rather than pooled naively across all rows,
which would conflate real by-design heterogeneity across cells with
genuine Monte Carlo error; see
`results/benchmarks/confirmation/compute_confirmation_mcse.R` for the
reproducing code and `study1_confirmation_mcse.csv` for every cell - is
under 0.75 percentage points on three of the four primary metrics below
(Type I error, spurious-edge rate, false augmentation-only significance
rate) and 1.01 percentage points on the fourth and widest, power
recovery, reaching this precision at a small fraction of the compute a
full-grid, 100-replication confirmation would have cost, without changing
any conclusion reported here. The original pilot is retained as
exploratory context in `results/benchmarks/FINDINGS.md` in the package
repository; the confirmation run's full results and reproducing scripts
are reported in `results/benchmarks/confirmation/` there. Because
`sydni_brm`'s confirmed grid is narrower than the other four SYDNI
engines', any direct comparison against it below uses the matching subset
of their data rather than their full grid, and is noted as such - the
naive full-grid comparison mechanically makes `sydni_brm` look worse on
power recovery than it actually is, an artifact of scenario/sample-size
mix rather than a real difference.

**Off-label stress test.** Every SYDNI engine and `synthpop`'s CART
variant clearly outperformed the off-label bootstrap/SMOTE/ADASYN
comparators on Type I error, consistently across all three sample sizes
tested: a pooled Type I error rate of 0.162 (\(\pm\)0.001 95% MCI
half-width, design-stratified by method/scenario/n\(_{\mathrm{obs}}\)/
augmentation-ratio cell - see
`results/benchmarks/confirmation/study1_confirmation_mcse.csv`) averaged
across the five SYDNI engines and `synthpop_cart`, against 0.248
(\(\pm\)0.002) averaged across bootstrap/noise-bootstrap/SMOTE/ADASYN -
roughly 1.5x, not the 2x suggested by the smaller pilot. This was true
even though these
comparator methods were evaluated under the same discouraged
pooled-analysis usage as SYDNI itself in this stress test, not under
their own intended use cases.

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
baseline rather than a real comparator (see below) -
but both remain in the underlying pooled results summarized in this
section.

**Intended-use safety.** Under SYDNI's documented secondary-check
workflow, the false-positive rate for effects that were null in the
observed data but spuriously significant only after augmentation
(`false_augmentation_only_sig_rate`) was lowest among the four
continuous-outcome SYDNI engines for `sydni_graphical` (0.102,
\(\pm\)0.002 95% MCI half-width), and highest for `sydni_copula` (0.132,
\(\pm\)0.003, narrowly above `sydni_bart`'s 0.126, \(\pm\)0.003) - all
still well below every off-label comparator (0.167-0.204). (`sydni_factor`'s
rate on this metric, computed on its discretized outcome and not directly
comparable to the figures above, is reported separately below: 0.087 -
see the `sydni_factor` subsection.)

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

**Structure retention and engine tradeoffs.** A raw mean absolute
correlation-matrix error is deliberately not used as the headline
structure metric here (it is reported in Appendix B instead): it is
mechanically larger for any engine that deliberately shrinks correlations
toward a conservative target as part of its false-positive protection, so
an engine doing exactly what it is designed to do looks unconditionally
"worse" on that number than one that does not shrink at all - a higher
value is not evidence of a defect, and can be the opposite. Two metrics
below separate that conflation: `null_edge_spurious_rate`, the share of
pairs the observed data already showed no meaningful relationship for
(\(|r| < 0.10\)) that the augmented data manufactured into an apparent
correlation (\(|r| > 0.20\)) - a direct test of whether an engine's noise
injection avoids inventing spurious structure - and
`cor_sign_reversal_rate`, restricted to pairs with a real observed
relationship, which asks whether that relationship's *direction* survived
augmentation. No single SYDNI engine dominated on every axis; the five
engines were designed as different tradeoffs rather than competitors on
one scale, and the confirmation run's role is to show what each tradeoff
actually costs and buys.

Because `sydni_brm` was confirmed on a narrower grid (above), the other
engines below are restricted to the matching subset of scenarios, sample
sizes, and augmentation ratios rather than their full confirmed grid, so
every row in this table is drawn from a like-for-like comparison.
`sydni_factor` is deliberately excluded from this table and reported
separately afterward: `engine_factor()` requires ordinal/integer input,
so its numbers come from a discretized version of the continuous outcome
`y`, while every engine below uses the same continuous `y` used
elsewhere in this benchmark - ranking it alongside engines evaluated on
a different outcome would imply a like-for-like comparison that does not
hold.

\begin{center}
\begin{tabular}{lrrrr}
\hline
Engine & Type I error (rank) & Power recovery (rank) & Spurious-edge rate (rank) & sec/condition \\
\hline
\texttt{sydni\_graphical} & 0.143 (1st, safest) & 0.053 (4th, worst)   & 0.0008 (1st, best)  & 0.03 \\
\texttt{sydni\_brm}       & 0.148 (2nd)         & 0.059 (3rd)          & 0.0009 (2nd)        & 113.21 \\
\texttt{sydni\_copula}    & 0.199 (3rd)         & 0.080 (1st, best)    & 0.0019 (3rd)        & 0.03 \\
\texttt{sydni\_bart}      & 0.203 (4th, worst)  & 0.073 (2nd)          & 0.0150 (4th, worst) & 0.30 \\
\hline
\end{tabular}
\end{center}

At this replication count (1,000-2,000 conditions per engine on the
matched grid), most of these ranks are resolved well outside Monte Carlo
noise - the design-stratified 95% MCI half-width (full per-cell detail in
`results/benchmarks/confirmation/study1_confirmation_mcse.csv`) is under
0.8 percentage points on every cell in this table except `sydni_brm`'s
power recovery (\(\pm\)1.01 percentage points, this table's widest cell).
Three adjacent pairs are close enough that their rank order should be
read as directional rather than confirmed, each roughly one combined
standard error apart: `sydni_graphical` and `sydni_brm` on Type I error
(0.143 vs.\ 0.148) and on power recovery (0.053 vs.\ 0.059); and
`sydni_copula` and `sydni_bart` on Type I error (0.199 vs.\ 0.203).
`sydni_bart`'s spurious-edge rate is the one large, unambiguous outlier
in this table - roughly 8-19x every other engine's rate, including
`sydni_brm`'s, and resolved well outside any of these margins - and is
discussed further below.

The runtime column changes the interpretation of this table substantially.
`sydni_brm` is never worst on any of the three axes above (2nd, 3rd, 2nd)
- a real, defensible sense of "balanced" - but costs roughly three orders
of magnitude more compute than every other engine. Comparing each engine
against `sydni_brm` directly (gap = engine minus `sydni_brm`; negative
Type I error/spurious-rate gap and positive power gap both indicate the
cheaper engine performing better) shows how much of that balance survives
at a fraction of the cost:

\begin{center}
\begin{tabular}{lrrrr}
\hline
Engine & Speedup vs.\ \texttt{sydni\_brm} & Type I gap & Power gap & Spurious-rate gap \\
\hline
\texttt{sydni\_graphical} & 3{,}789x & $-$0.005 (\textasciitilde tied) & $-$0.006 (worse) & $-$0.0001 (\textasciitilde tied) \\
\texttt{sydni\_copula}    & 3{,}580x & $+$0.051 (worse)  & $+$0.022 (better)               & $+$0.0010 (worse) \\
\texttt{sydni\_bart}      & 380x     & $+$0.056 (worse)  & $+$0.014 (better)  & $+$0.0141 (much worse) \\
\hline
\end{tabular}
\end{center}

This splits the (non-`sydni_factor`) family into three roles rather than
a single ranking. `sydni_brm` is the reference point when compute is not
a binding constraint and no single axis dominates the analyst's
priority, at a cost of the rare small/sparse-data ADVI fitting failure
noted below. `sydni_graphical` is the strongest all-around cheap
alternative on the fair structure metric: it leads these four engines on
both Type I error and spurious-edge suppression at 3,789x less compute
than `sydni_brm`, giving up ground only on power recovery, where it is
the family's weakest. `sydni_copula` remains the power-recovery leader by
a clear margin and is the better choice when recovering real effects
matters more than minimizing false positives, at the price of a Type I
error gap comparable to `sydni_bart`'s. `sydni_bart` is the weakest
engine in the family on the fair structure metric - its spurious-edge
rate is 8-19x every other engine's, including `sydni_brm`'s - and is not
recommended as a default choice; its one genuine advantage is power
recovery competitive with `sydni_copula`'s at less compute, which may
still make it worth reaching for when that specific tradeoff is wanted
and the elevated spurious-correlation rate is an acceptable cost. No
engine needs to lead every metric to be useful: no simulation-based
benchmark can realistically cover the full space of real-world data an
analyst might encounter, and having several structurally different,
inexpensive options available is itself part of the value of SYDNI's
multi-engine design, independent of any one engine's ranking here.

**`sydni_factor`'s results: ordinal outcome, reported separately.**
Because `engine_factor()` was evaluated against a discretized version of
`y` rather than the continuous outcome every other engine above used,
its numbers are not placed in the ranked table, the speedup table, or
any "Nth of five" comparison - doing so would imply a like-for-like
comparison that the design of this benchmark does not support. On the
same matched grid used above, `sydni_factor`'s own numbers (95% MCI
half-width in parentheses; see `study1_confirmation_mcse.csv`) are: Type
I error 0.129 (\(\pm\)0.005), power recovery 0.062 (\(\pm\)0.007),
spurious-edge rate 0.0006 (\(\pm\)0.0002), at 1.28 seconds per condition;
on the full (non-matched) grid, its false augmentation-only significance
rate - the intended-use safety metric discussed above - is 0.087
(\(\pm\)0.002). Within its own ordinal-outcome frame, these numbers
are directionally favorable - lower Type I error and a lower spurious-edge
rate than any of the four continuous-outcome engines above, at a runtime
in the same broad cheap-engine range (fractions of a second to a few
seconds) - but whether that favorable read would hold under the same
continuous outcome the rest of this study uses is an open question this
benchmark cannot answer, not a confirmed finding. Analysts whose actual
outcome is genuinely ordinal or discretized should read this row as the
directly relevant evidence; analysts comparing engines for a continuous
outcome should use the four-engine table above instead.

Separately from the cross-engine comparison above, `sydni_brm` showed
zero fitting failures across all 1,000 confirmed conditions in this
study - an improvement over the original 10-replication pilot's 4
failures in 630 conditions, though the underlying ADVI fragility that
caused those pilot failures is real and was confirmed to recur, at a
similar rate, in the Section 8.7 measurement-aware confirmation run (2 of
400 conditions; see `KNOWN_ISSUES.md` in the package repository) - at a
computational cost roughly three orders of magnitude larger than any
other engine (approximately 113 seconds per condition versus fractions
of a second).

![Figure 3: Mean runtime per condition by method, log scale.](figures/figure3_diagnostics.png)

**Figure 3.** Mean runtime per condition by method (log10 scale),
SYDNI engines against off-label comparators, confirmation-run data.
`observed_only` and `synthpop`'s sample-based variant are omitted as
near-zero-runtime reference/floor baselines rather than real comparators
(see above). One degenerate simulated draw - `small_cell`, \(n_{\mathrm{obs}} = 75\),
\(n_{\mathrm{syn}}/n_{\mathrm{obs}} = 2\), replicate 97, where the
smallest-cell subgroup happened to receive zero observations, collapsing
`group` to a constant column - affected 11 of the 69,300 non-`sydni_brm`
method-level conditions (one row per method that fit on that replicate),
now handled as 11 recorded failures rather than a benchmark-harness
crash - see `KNOWN_ISSUES.md`; 0 of 1,000 `sydni_brm` conditions failed
(see above).
Every method also emits occasional informational package notices (e.g.,
ordinal-as-numeric, small-subgroup, or missing-`measurement_blocks`
notices) unrelated to
fitting problems; the raw counts are reported alongside the reproducing
scripts in `results/benchmarks/FINDINGS.md`.

**Sample-size dependence.** Two patterns emerged from varying
\(n_{\mathrm{obs}}\), one confirmed and one revised at this replication
count. First, the power-recovery rate declined toward zero for every
method as \(n_{\mathrm{obs}}\) increased (Figure 2), consistent with
augmentation mattering primarily when the observed sample is itself
underpowered - a pattern that supports rather than undermines SYDNI's
intended positioning, and holds cleanly at 100 replications. Second, the
pilot's claim that Type I error for null coefficients increased with
\(n_{\mathrm{obs}}\) "for every method, including the observed-only
baseline" does not survive confirmation: `observed_only` is essentially
flat across \(n_{\mathrm{obs}}\) (0.064/0.067/0.064), and the pattern
splits cleanly by family rather than being universal - it holds for
`sydni_factor`, `sydni_graphical`, `sydni_copula`, and `synthpop_cart`,
but the four off-label comparators and `sydni_bart` trend flat-to-*
decreasing* instead (Figure 1). Since `observed_only` never touches
augmentation, whatever drives the increasing group's pattern is a
property of those specific augmentation methods interacting with sample
size, not a generic analysis-model-misspecification effect across every
method as the pilot's smaller sample suggested. This does not change the
relative method rankings, which were stable across sample sizes, but the
mechanism claim itself is revised rather than confirmed.

**Does BART show a distinct nonlinear/interaction advantage? A direct
test found none.** `sydni_bart`'s practical role - the weakest engine in
the family on the fair structure metric, with a power-recovery edge as
its one clear advantage - is addressed above. A
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

### 8.7 Measurement-aware calibration results

Section 8.6 tests inference on manifest predictor variables directly - the
variables synthesized are the variables analyzed. It never exercises the
measurement-aware calibration mechanism of Section 4.1
(`measurement_blocks` / `between_block_shrink`), which is specifically
designed to protect within-scale item reliability while shrinking
between-construct associations. A second, companion benchmark was run to
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
for reference. As with Section 8.6, an initial 10-replication pilot
motivated a confirmation run at higher replication, reported as the
primary result throughout this section: 100 replications per condition
for every method except `sydni_brm_flat`/`sydni_brm_ma` (13,500
conditions, the same full \(n_{\mathrm{obs}} \in \{75, 150, 300\}\) and
augmentation-ratio grid as Section 8.1), and 50 replications on a
narrower predeclared grid for the two `sydni_brm` variants
(\(n_{\mathrm{obs}} \in \{75, 150\}\), \(n_{\mathrm{syn}}/n_{\mathrm{obs}}
\in \{1, 2\}\), 400 conditions), for the same compute-cost reasons as
Section 8.6. 4 of these 400 `sydni_brm` conditions failed (all
`sydni_brm_ma`; 2 the same ADVI fitting crash documented in Section 8.6,
2 a transient compilation-time error confirmed non-reproducible on
isolated retry - see `KNOWN_ISSUES.md` in the package repository); every
other condition across both studies completed. The original pilot is
retained as exploratory context in
`results/benchmarks/simulation_study_2/FINDINGS.md` in the package
repository; the confirmation run's full results and reproducing scripts
are reported in `results/benchmarks/confirmation/` there. This section
summarizes the headline findings, organized around three scenarios this
design speaks to.

**Scenario 1: manifest/composite predictors with no indicator-level data.**
Not tested by this study - this is Section 8.6's scenario, and nothing
here revises its conclusions. The DGP and analysis target differ enough
(no latent/indicator/aggregation layer at all) that this study's numbers
do not generalize to that case in either direction.

**Scenario 2: indicator-level data exists, but only composite-score output
is needed.** Here the flat-vs-ma comparison is directly informative, and
the result holds across every engine, not just one: every `flat`
variant's Type I error and false discovery rate (0.188-0.248 /
0.105-0.138) sits in the same range as the off-label comparators
(0.217-0.245 / 0.117-0.130), while every `ma` variant (0.060-0.084 /
0.036-0.049) sits well below both, though not as close to the
observed-only baseline (0.051 / 0.032) as the smaller pilot suggested -
`ma`'s Type I error runs roughly 1.2-1.6x the baseline's rather than
nearly matching it. The practical implication is unchanged:
`measurement_blocks` is worth supplying whenever indicator-level data
exists, even when the end goal is only a composite-score regression -
skipping it costs real inferential safety on every engine tested, at a
runtime premium that stays under a quarter-second per condition for
`sydni_copula`, `sydni_graphical`, and `sydni_bart`, but rises to roughly
21x for `sydni_brm` (17.5 to 361.8 seconds).

![Figure 4: Flat vs. measurement-aware, paired by engine, across n_obs.](figures/figure4_flat_vs_ma.png)

**Figure 4.** Type I error rate, false discovery rate, and power recovery
rate for each of the five SYDNI engines run without `measurement_blocks`
(`flat`, left) against the same five engines run with the full
measurement-aware calibration path (`ma`, right), paired by engine and
plotted across \(n_{\mathrm{obs}}\), confirmation-run data. `sydni_brm`'s
line ends at \(n_{\mathrm{obs}} = 150\) in both panels, reflecting its
narrower confirmed grid (see above). Panel/label/reference-band
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
\texttt{sydni\_factor\_ma}    & 2.00 & 0.059 & 0.035 & 0.048 & 0.037 & 0.031 & 0.52 \\
\texttt{sydni\_graphical\_ma} & 2.17 & 0.061 & 0.036 & 0.058 & 0.041 & 0.036 & 0.22 \\
\texttt{sydni\_copula\_ma}    & 3.00 & 0.064 & 0.038 & 0.051 & 0.063 & 0.032 & 0.23 \\
\texttt{sydni\_bart\_ma}      & 3.33 & 0.074 & 0.043 & 0.083 & 0.100 & 0.029 & 0.74 \\
\texttt{sydni\_brm\_ma}       & 4.50 & 0.074 & 0.043 & 0.046 & 0.049 & 0.055 & 361.25 \\
\hline
\end{tabular}
\end{center}

**All four non-`sydni_brm` rows above are computed on the same
\(n_{\mathrm{obs}} \in \{75, 150\}\), \(n_{\mathrm{syn}}/n_{\mathrm{obs}}
\in \{1, 2\}\) subset that `sydni_brm_ma` was actually confirmed on** -
not the full three-point grid used elsewhere in this section - so that
`sydni_brm_ma` is compared on equal footing rather than against numbers
drawn from cells it was never run on. An earlier draft of this table
mixed the two grids, which understated `sydni_bart_ma`'s power recovery
(0.052 vs.\ 0.083 here) and, more importantly, made `sydni_brm_ma` look
more competitive than the matched comparison supports.

**Precision caveat: only 4 design cells support this table.** Unlike
Section 8.6's 20-cell matched grid, this scenario's grid is
\(n_{\mathrm{obs}} \in \{75, 150\} \times n_{\mathrm{syn}}/n_{\mathrm{obs}}
\in \{1, 2\}\) - 4 design cells - so the design-stratified 95% MCI
half-width on the rate-based columns (Type I error, FDR, power recovery,
between-construct sign-reversal) is substantially wider than Section
8.6's: roughly 1.0-2.6 percentage points depending on engine and metric,
versus under 0.8pp there (full per-cell detail in
`results/benchmarks/confirmation/study2_confirmation_mcse.csv`,
reproduced by `compute_confirmation_mcse.R` in the same directory). At
that precision, essentially none of the pairwise differences among the
five engines on Type I error, FDR, or power recovery are resolved outside
roughly one combined standard error - the mean-rank ordering below on
those three axes should be read as suggestive, not statistically
confirmed. Within-scale correlation error and runtime are the exception:
their half-widths are small relative to the between-engine gaps (e.g.\
`sydni_bart_ma`'s within-scale error, 0.100 \(\pm\)0.001, is clearly
separated from every other engine's 0.037-0.063, half-widths
0.0004-0.0055; `sydni_brm_ma`'s runtime, 361 \(\pm\)14 seconds, is
separated from every other engine's by orders of magnitude) - those two
conclusions rest on firmer footing than the rate-based ranking.

Subject to that caveat, one finding here revises this family's Section
8.6-derived positioning, scoped strictly to Scenario 3. First, by point
estimate, **`sydni_brm`'s measurement-aware path is not the all-arounder
its manifest-variable behavior would suggest.** Where Section 8.6 found
`sydni_brm` never worst on any axis, here it finishes last by mean rank
and is worst on four of the six axes - Type I error, power recovery,
between-construct sign-reversal rate, and runtime (the last at roughly
470-1,600x the compute of the other four engines, resolved well outside
Monte Carlo noise) - while sitting mid-pack on false discovery rate (4th)
and within-scale correlation error (3rd, also resolved). Of these,
runtime and within-scale error are confirmed weaknesses; Type I error,
power recovery, and between-construct sign-reversal rate are directional
at this replication count, per the precision caveat above. This is a
broader weakness than the
grid-mismatched comparison in an earlier draft suggested (that version
reported `sydni_brm_ma` worst on only two of six axes and second-best on
power recovery); on the matched grid it is worst on power recovery too.
Notably, the smaller 10-replication pilot's original characterization -
worst or tied-worst on four of six axes - turns out to be closer to this
matched-grid picture than the grid-mismatched full comparison was. This
does not touch Section 8.6's manifest-variable conclusion about
`sydni_brm`; it is confined to the measurement-aware use case. Second, by
point estimate, **`sydni_bart`'s niche narrows for this scenario without
being eliminated.** It still posts the best power recovery (0.083) and the best
between-construct sign-reversal rate (0.029) of the five, but its
within-scale correlation error (0.100) is roughly 1.6-2.7x every other
`ma` engine's (0.037-0.063) - the one engine whose measurement-aware path
does not clearly protect the within-construct reliability the mechanism
is specifically designed to protect. Consistent with Section 8.6's own
finding that `sydni_bart` is not recommended as a default choice, for
measurement-aware use it is better read as a specialist for when power or
between-construct fidelity matters more than item-level reliability, not
as a default. `sydni_factor_ma` and `sydni_graphical_ma` remain the
recommended defaults for this scenario, though on the matched grid they
are no longer tied and neither leads on every axis: by point estimate,
`sydni_factor_ma` posts the best mean rank overall (2.00) and is best or tied-best on Type
I error, false discovery rate, and within-scale correlation error, at the
cost of the weakest power recovery among the four non-`sydni_brm` engines
(4th of 5 overall). `sydni_graphical_ma` (mean rank 2.17) is faster still
(0.22 vs.\ 0.52 seconds per condition) and recovers more power (2nd of 5),
at the cost of the weakest between-construct sign-reversal rate among the
four non-`sydni_brm` engines (4th of 5 overall). The choice between them
is a safety-leaning default (`sydni_factor_ma`) versus a faster,
power-leaning one (`sydni_graphical_ma`), not a clear dominance of one
over the other.

Combining all three scenarios gives a scenario-dependent recommendation
rather than one global engine ranking:

\begin{center}
\begin{tabular}{p{4.2cm}p{5.3cm}p{2.8cm}}
\hline
User need & Recommended engine(s) & Evidence source \\
\hline
Manifest/composite predictors, no indicator data & \texttt{sydni\_brm} (reference); \texttt{sydni\_graphical} (strongest cheap all-arounder); \texttt{sydni\_copula} (power-recovery leader); \texttt{sydni\_bart} (power specialist, not a default); \texttt{sydni\_factor} (ordinal outcomes only - not cross-engine comparable) & Section 8.6 \\
Indicator-level data, composite-score output only & Supply \texttt{measurement\_blocks} regardless of engine & Section 8.7 (flat vs.\ ma) \\
Indicator-level data, item-level/SEM output needed & \texttt{sydni\_factor\_ma}/\texttt{sydni\_graphical\_ma} (defaults); \texttt{sydni\_bart\_ma} (power/between-construct specialist); \texttt{sydni\_brm\_ma} (viable but expensive, no demonstrated payoff) & Section 8.7 (ma arm) \\
\hline
\end{tabular}
\end{center}

**Sample-size dependence** mirrors Section 8.6's pattern on this different
DGP: power recovery for every non-`sydni_brm` `ma` engine shrinks toward
zero as \(n_{\mathrm{obs}}\) grows (e.g. `sydni_factor_ma`: 0.048 to
0.005 from \(n_{\mathrm{obs}} = 75\) to 300; `sydni_brm_ma` shows the same
direction over its narrower confirmed range, 0.077 to 0.015 from
\(n_{\mathrm{obs}} = 75\) to 150), and each `ma` engine's Type I error
narrows toward, without fully reaching, the observed-only baseline's own
rate by \(n_{\mathrm{obs}} = 300\) (baseline 0.053; `ma` engines
0.062-0.097, still 1.2-1.8x the baseline), consistent with SYDNI's
positioning as a small/underpowered-sample tool rather than a
general-purpose one.

**Open questions carried forward from this study.** Why `sydni_bart`'s
measurement-aware path does not protect within-scale reliability the way
the other four engines' do is not diagnosed here - a documented open
question alongside Section 8.6's still-open BART nonlinearity question,
not a resolved one. `sydni_brm_ma`'s runtime cost was traced during this
study's setup to a `brms`/Stan mechanics issue rather than a modeling
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
the workflow be illustrated with a real psychometric instrument rather than
only a designed data-generating process. This section reports a single-run
illustration (Study 3), not a
replicated benchmark: real data has no known ground truth, so the metrics
below are structure-retention and downstream-model-stability comparisons
rather than the Type I error/power-recovery metrics Sections 8.6-8.7 use
against known truth. The full 17-method exploratory run, its raw outputs, and
the additional engine comparisons are retained in the supplementary results
document, `paper/sydni-jss-supplementary-results.md`.

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

The primary sensitivity set uses the measurement-aware `copula`, `factor`, and
`graphical` engines. These choices reflect the ordinal item-battery structure:
they respectively represent dense joint dependence, latent measurement
structure, and sparse direct dependence. Each run uses `n_syn = n_obs` and the
same declared measurement blocks. The full five-engine, flat-versus-
measurement-aware, and off-label-comparator runs are supplementary rather than
part of this worked analyst workflow.

**Downstream-model stability.** The observed model,
\(\mathrm{Grit} \sim E + N + A + O\), finds Neuroticism (b = -0.30,
p < .001) and Agreeableness (b = 0.16, p < .001) significant, with
Extraversion (b = 0.005, p = .88) and Openness (b = 0.076, p = .10) not
significant. The measurement-aware copula and graphical runs preserve that
pattern. The factor run leaves the Openness point estimate essentially
unchanged (0.076 to 0.076) but yields p = .038. Under the reporting protocol
in Section 7, this synthetic-only significance change is exploratory: it
identifies sensitivity to the synthesis choice rather than additional evidence
for an Openness association.

**Structure diagnostics.** The three primary runs preserve within-scale
structure closely and all reduce between-construct correlation magnitude, as
requested by `between_block_shrink = 0.5`.

| Engine | Within-scale alpha change | Within-scale correlation error | Between-construct magnitude change |
|---|---:|---:|---:|
| `copula` | 0.016 | 0.030 | -0.054 |
| `graphical` | 0.018 | 0.032 | -0.054 |
| `factor` | 0.011 | 0.036 | -0.055 |

These values describe one generated extension, not a general comparative
ranking. They demonstrate the diagnostics an analyst should inspect before
interpreting a downstream sensitivity result.

This illustration is a worked example, not an independent confirmation of
Sections 8.6-8.7. A single real dataset with no replication cannot carry the
same evidential weight as a simulation study with known truth. Disagreement
across the selected engines is useful here because it makes explicit which
conclusions depend on a particular plausible generative assumption.

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
methodology vignettes under `vignettes/`. Before submission, the article must
be accompanied by a versioned replication bundle that regenerates every table
and figure from raw inputs, records `sessionInfo()` and package versions, and
documents runtime and optional-engine requirements. The submitted manuscript
must likewise be rendered in the Journal of Statistical Software LaTeX style;
this Markdown file is the working source rather than the submission artifact.

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
   model assumptions. The pilots did not detect a BART advantage under their
   specific nonlinear design, and the measurement-aware pilot showed a
   trade-off between BART's structure retention and downstream metrics. These
   results are design-specific rather than general engine properties.
8. The computational cost differs substantially across engines, especially for
   BRM models requiring Stan backends. Runtime should be measured in the
   analyst's own setting rather than inferred from the pilot alone.
9. In BRM's measurement-aware mode, structural-family-score calibration is
   applied before item-level donor reconstruction, so final item-level
   composite correlations can differ slightly from the calibration target,
   unlike the other four engines.
10. Retention-aware group balancing currently supports only
    `balance_target = "largest_group"`.

These limitations are consistent with the intended use of SYDNI as a
robustness workflow, not a universal synthetic-data generator.

## 12. Discussion

SYDNI frames synthetic augmentation as a secondary sensitivity-analysis tool.
The key question is not whether synthetic rows should be treated as additional
observations, but whether analytic conclusions are stable when the observed
sample is extended through controlled perturbation and structure-preserving
synthesis. The original-data analysis remains primary. Analysts should
pre-specify a defensible set of engines, settings, and seeds; compare direction,
effect magnitude, uncertainty, and substantive decisions across those runs; and
treat synthetic-only findings as exploratory signals rather than confirmation.

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

The two benchmark studies provide targeted validation evidence, not a global
ranking of engines. Both are now confirmed at substantially higher
replication than their originating pilots (100 per condition for every
method except `sydni_brm`/`sydni_brm_flat`/`sydni_brm_ma`, 50 on a
predeclared focused grid for those, in each study). The manifest-variable
study found a robust separation from off-label resampling comparators and
did not detect a BART advantage in its nonlinear design. The
measurement-aware study evaluated the calibration mechanism under one
item-battery design and found `measurement_blocks` worth supplying
whenever indicator-level data exists. Together, these results inform
engine-specific cautions and motivate further testing; they do not
establish a universal winner or a formal inferential guarantee.

A single-run illustration on the package's `grit` dataset (Section 9) shows
how the workflow can be reported with item-level data. Its structure-retention
and downstream-model summaries are descriptive: the example has neither known
truth nor replication, and it cannot independently validate the benchmark
results.

Future work should address external validity now that both studies carry
adequate Monte Carlo precision on their own designs: additional
data-generating designs and empirical examples are needed, for both
studies, to test generalizability beyond the specific scenarios and
item-battery structure exercised here. Other priorities include improved
synthetic-data diagnostics, missing-data integration, and optional
disclosure-risk diagnostics
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

## Appendix B: Raw correlation-matrix error, by method

Section 8.6 deliberately does not use raw mean absolute correlation-matrix
error (`mean_abs_cor_error`) as a headline structure-retention metric: it
is mechanically larger for any method that deliberately shrinks
correlations toward a conservative target, so a method doing exactly what
it is designed to do looks unconditionally "worse" on this number than
one that does not shrink at all. It is reported here for completeness,
confirmation-run data (100 replications for every method except
`sydni_brm`, full grid; 50 replications for `sydni_brm` on its narrower
confirmed grid - Section 8.6):

\begin{center}
\begin{tabular}{lr}
\hline
Method & Mean absolute correlation error \\
\hline
\texttt{bootstrap}        & 0.028 \\
\texttt{ADASYN}           & 0.027 \\
\texttt{noise\_bootstrap} & 0.031 \\
\texttt{SMOTE}            & 0.032 \\
\texttt{sydni\_copula}    & 0.038 \\
\texttt{synthpop\_cart}   & 0.044 \\
\texttt{sydni\_bart}      & 0.051 \\
\texttt{sydni\_brm}       & 0.065 \\
\texttt{sydni\_graphical} & 0.071 \\
\texttt{sydni\_factor}    & 0.076 \\
\hline
\end{tabular}
\end{center}

Read against Section 8.6's fair structure metrics
(`null_edge_spurious_rate`, `cor_sign_reversal_rate`), the ordering here
is close to inverted for several methods: `sydni_graphical` posts the
second-largest raw correlation error in this table while simultaneously
posting the lowest spurious-edge rate of the four continuous-outcome
engines in Section 8.6's engine-tradeoff table (0.0008); `sydni_factor`
posts the single largest raw correlation error here while posting the
lowest spurious-edge rate on its own ordinal-outcome terms in the
`sydni_factor` subsection of Section 8.6 (0.0006, not directly comparable
to `sydni_graphical`'s figure - see that subsection); and the four
off-label comparators post the smallest raw errors here while showing no
deliberate shrinkage mechanism at all - they simply perturb the data
little enough that little moves. This table should be read as a
description of how much each method's synthetic data numerically differs
from the observed correlation structure, not as a ranking of data
quality.

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
