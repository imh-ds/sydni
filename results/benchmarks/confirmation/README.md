# Confirmation results: reproduction from saved simulations

This folder contains per-replication summaries from two synthetic simulation studies and the scripts needed to reconstruct the retained Monte Carlo summaries and paired-block contrasts. The CSV inputs are simulation outputs; they are not survey records or additional human observations. Four raw files correspond to study 1 and study 2, each split into Bayesian regression (BRM) and other methods.

## What you can reproduce

- `compute_confirmation_mcse.R` reads all four raw CSVs and writes the two `study*_confirmation_mcse.csv` summaries. These are retained historical design-stratified calculations. The original calculations treated design cells as independent despite shared simulation streams and are not the paired intervals used in the methodology evaluation. Do not use them as precision-based rankings.
- `02_paired_block_reanalysis.py` reads three raw CSVs and reconstructs `paired_block_bootstrap.csv` using complete replicate blocks, 5,000 bootstrap samples, and seed 20261105. The paired contrasts are conditional on the specified grids and complete blocks.

The local confirmation runs recorded R 4.6.0 on Windows and package commit `d64fd74a9d877eb7626d226046e09eb82d8c5616`. The original benchmark frameworks were untracked at that commit, so this repository does not establish their exact historical executed contents. These instructions reconstruct saved-result summaries; they do not recreate the original simulation runs or recover their complete environment. The BRM runs used a smaller grid and had recorded failures.

## Run in a disposable checkout

From the repository root, with R 4.6.0 and Python 3.14.3 available, run:

```powershell
Rscript results/benchmarks/confirmation/compute_confirmation_mcse.R
$env:SYDNI_CONFIRMATION_DIR = 'results/benchmarks/confirmation'
$env:SYDNI_AUDIT_OUT = 'confirmation-reproduced'
python results/benchmarks/confirmation/02_paired_block_reanalysis.py
```

The R script uses base R and overwrites the two MCSE CSVs in this folder; use a disposable checkout or copies. The Python script uses only the standard library and writes `confirmation-reproduced/paired_block_bootstrap.csv`; compare it with the saved `paired_block_bootstrap.csv` here. In a local check, both R outputs and the Python output matched all corresponding reference CSV rows and fields exactly. No dependencies were installed or simulations refitted for that check. The `SHA256SUMS.txt` manifest records the repository blob bytes. Git preserves the raw input CSVs and paired output without line-ending conversion; other text files can differ in a checkout with automatic line-ending conversion. File hashes do not establish historical provenance.
