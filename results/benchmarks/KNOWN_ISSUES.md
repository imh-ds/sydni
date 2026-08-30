# Known issues

## small_cell: one degenerate draw during the Study 1 confirmation run

11 of 69,300 rows in `confirmation/study1_nonbrm_confirmation_raw.csv`
(~0.02%) have `failed == TRUE`, all from a single simulated draw:
`scenario = small_cell`, `n_obs = 75`, `aug_ratio = 2`, `replicate = 97`
- one row per method that attempted to fit on that replicate (`sim_seed`
depends on `replicate`/`n_obs`/`aug_ratio` but not `method`, so every
method independently drew the same degenerate data).

Root cause: at `n_obs = 75`, `small_cell`'s smallest subgroup can by
chance receive zero observations, collapsing the `group` covariate to a
constant column. `lm()`'s contrasts then throw `"contrasts can be
applied only to factors with 2 or more levels"` (or an equivalent
degenerate-fit error depending on the method). This is the same failure
mode described in the `error_row()` comment at
`inst/benchmarks/inferential_safety_framework.R:608-618`, which explains
why every method's fit on this replicate is caught and recorded as a
graceful failure instead of crashing the parallel worker.

Not treated as a bug to fix: it's an expected, low-probability
consequence of drawing a subgroup indicator at small `n_obs`, not a
package defect. It is caught and recorded as 11 graceful `failed = TRUE`
rows (via the `error_row()`/backstop-`tryCatch` pattern) rather than
crashing the run, so it costs 11 lost method-level conditions rather
than corrupting or blocking anything else. Left as-is; not blocking
interpretation of the remaining 69,289 rows.

## sydni_brm: rare ADVI fitting failures on small/sparse data

4 of 630 `sydni_brm` rows in `inferential_safety_pilot_raw.csv` (~0.6%) have
`failed == TRUE`, all at `n_obs = 75`, `aug_ratio = 2`: 2 of 4 in the
`small_cell` scenario (replicates 4 and 9), and 1 each in `null_only`
(replicate 4) and `moderation_imbalance` (replicate 4).

Reproduced deterministically (same seed -> same failure), so this is a real
numerical failure mode, not run-to-run flakiness. The underlying error:

```
Warning: Fitting finished unexpectedly! Use the $output() method for more information.
...
Error: Fitting failed. Unable to retrieve the metadata.
```

cmdstanr's fullrank-ADVI optimizer crashes outright (not a convergence
warning) in the smallest-sample, largest-augmentation-ratio corner of this
grid. The failures span three scenarios, so this pilot does not isolate
sparsity as their cause. Consistent with brms/`algorithm = "fullrank"`'s known reliability
limits on small or sparse designs (see the earlier Pareto-k reliability
discussion for `sydni_brm`).

Not treated as a bug to fix: this matches known, expected guidance for small
subgroups - train on the larger sample rather than pre-filtering to the
exact small subgroup of interest before fitting. Left as-is; not blocking
interpretation of the rest of the `sydni_brm` results, which are otherwise
0% failures across the remaining 626 rows.

## sydni_brm_ma: rare failures during the Study 2 confirmation run

4 of 400 `sydni_brm_ma`/`sydni_brm_flat` rows in
`confirmation/study2_brm_confirmation_raw.csv` (~1%) have `failed == TRUE`,
all `sydni_brm_ma` (none `sydni_brm_flat`), of two distinct kinds:

**2 of 4: the same ADVI fitting crash documented above**
(`n_obs = 150, aug_ratio = 2, replicate = 20` and
`n_obs = 75, aug_ratio = 1, replicate = 39`), same
`Error: Fitting failed. Unable to retrieve the metadata.` signature. Not a
new failure mode - the same known fullrank-ADVI reliability limit,
recurring here because `sydni_brm_ma` fits multiple per-construct BRM
submodels per condition, giving more chances for any one of them to hit it.

**2 of 4: a transient compilation-time error**
(`n_obs = 75, aug_ratio = 1, replicate = 17` and
`n_obs = 150, aug_ratio = 1, replicate = 26`), message
`An error occured during compilation! See the message above for more
information.` Investigated and confirmed **not reproducible**: both
conditions were re-run in isolation (single worker, no concurrent access to
the shared `stan_cache` directory) and succeeded cleanly on the first
retry, with no code or data change. The underlying model structure is
deliberately held fixed across replicates specifically to keep the
compiled Stan binary stable (see the `predictor_method = "full"` comment in
`simulate_measurement_aware_data()`), which rules out a data-dependent
degenerate-model explanation like the `small_cell` constant-covariate
draw documented above for Study 1 (recorded gracefully, not a bug to
fix). The most likely cause is
transient Windows-level contention (file locking or antivirus scanning)
when multiple parallel workers concurrently execute the same already-
compiled `.exe` - the same general class of Windows/cmdstanr fragility
already noted above for the persistent Stan model cache, just surfacing at
execute time instead of first-compile time. The compile warm-up pass added
to the confirmation runners prevents the first-compile race; it does not
guarantee every later concurrent execution is contention-free.

Not treated as a bug to fix: both failures are non-deterministic (unlike
the ADVI crash, which reproduces on the same seed every time) and are
already recorded gracefully as `failed = TRUE` rows with `error_message`/
`error_time` rather than crashing the run, so they cost one lost
replicate each rather than corrupting or blocking anything else. Left
as-is; not blocking interpretation of the remaining 396 rows.
