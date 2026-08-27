# Known issues

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
