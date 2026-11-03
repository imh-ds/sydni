# sydni 0.2.0

Summarized from the commit history; earlier work is not itemized.

## Engines and calibration

* Added the measurement-aware synthesis path: item blocks with
  `between_block_shrink`, available for the copula, graphical, BART, factor
  and BRM engines.
* Added retention-aware augmentation handling and structural score
  calibration for the BRM engine.
* Added a two-stage layer for the copula engine, and support for singleton
  copula shrinkage without measurement blocks.
* `brms` is now optional for the `brm` engine, and non-Gaussian BRM intercept
  priors are handled.
* Fixed copula noise seeding.

## Documentation and repository

* README gained a repository layout section and a license note; a short
  tutorial was added at `docs/tutorial.md`.
* Vignettes, benchmark scripts, `LICENSE` (GPL-2 text) and `CITATION.cff` are
  now tracked.
* Simulation confirmation studies are summarized under `results/benchmarks/`;
  known limitations of those studies are listed in
  `results/benchmarks/known_issues.md`.
