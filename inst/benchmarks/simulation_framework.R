benchmark_sydni_framework <- function(
    data,
    variables,
    covariates = NULL,
    n_aug = nrow(data),
    seed = 123
) {
  observed <- data[, c(variables, covariates), drop = FALSE]

  synth_copula <- sydni(
    data = data,
    variables = variables,
    covariates = covariates,
    n_aug = n_aug,
    engine = "copula",
    seed = seed
  )$data

  bootstrap_only <- data[sample(seq_len(nrow(data)), size = n_aug, replace = TRUE), c(variables, covariates), drop = FALSE]

  list(
    observed = observed,
    synth_copula = synth_copula,
    bootstrap_only = bootstrap_only
  )
}
