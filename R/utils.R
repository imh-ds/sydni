#' Internal utilities for SYDNI
#'
#' Helper utilities for reproducibility, input validation, and output
#' reconstruction.
#'
#' @keywords internal
#' @name sydni-internal
NULL

with_seed <- function(seed, code) {
  if (length(seed) == 1 && is.finite(seed) && !is.na(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      had_seed <- TRUE
    } else {
      old_seed <- NULL
      had_seed <- FALSE
    }

    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)

    set.seed(as.integer(seed))
  }

  force(code)
}

validate_engine_control <- function(engine_control) {
  if (is.null(engine_control)) {
    return(list())
  }

  if (!is.list(engine_control)) {
    stop("`engine_control` must be NULL or a named list.")
  }

  engine_control
}

resolve_sydni_columns <- function(data, variables = NULL, covariates = NULL) {
  data_names <- names(data)

  if (is.null(variables)) {
    variables <- setdiff(data_names, covariates)
  }

  if (length(variables) == 0) {
    stop("`variables` must select at least one column.")
  }

  if (!all(variables %in% data_names)) {
    missing_vars <- setdiff(variables, data_names)
    stop(
      paste0(
        "`variables` contains columns not present in `data`: ",
        paste(missing_vars, collapse = ", "),
        "."
      )
    )
  }

  if (!is.null(covariates) && !all(covariates %in% data_names)) {
    missing_covs <- setdiff(covariates, data_names)
    stop(
      paste0(
        "`covariates` contains columns not present in `data`: ",
        paste(missing_covs, collapse = ", "),
        "."
      )
    )
  }

  if (anyDuplicated(variables)) {
    stop("`variables` must not contain duplicates.")
  }

  if (!is.null(covariates) && anyDuplicated(covariates)) {
    stop("`covariates` must not contain duplicates.")
  }

  overlap <- intersect(variables, covariates)
  if (length(overlap) > 0) {
    stop(
      paste0(
        "Columns cannot be both synthesized and treated as covariates: ",
        paste(overlap, collapse = ", "),
        "."
      )
    )
  }

  list(
    variables = variables,
    covariates = covariates
  )
}

validate_sydni_data <- function(data, variables, covariates = NULL) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or tibble.")
  }

  if (nrow(data) == 0) {
    stop("`data` must contain at least one row.")
  }

  synth_data <- data[, variables, drop = FALSE]

  unsupported <- names(Filter(
    function(x) !(is.numeric(x) || is.integer(x) || is.factor(x) || is.ordered(x)),
    synth_data
  ))

  if (length(unsupported) > 0) {
    stop(
      paste0(
        "Unsupported synthesized column types: ",
        paste(unsupported, collapse = ", "),
        "."
      )
    )
  }

  zero_var <- names(Filter(
    function(x) length(unique(stats::na.omit(x))) < 2,
    synth_data
  ))

  if (length(zero_var) > 0) {
    stop(
      paste0(
        "Synthesized columns must have at least two distinct non-missing values: ",
        paste(zero_var, collapse = ", "),
        "."
      )
    )
  }

  if (!is.null(covariates)) {
    cov_data <- data[, covariates, drop = FALSE]
    bad_cov <- names(Filter(
      function(x) any(is.na(x)),
      cov_data
    ))

    if (length(bad_cov) > 0) {
      warning(
        paste0(
          "Covariates contain missing values and may create sparse groups: ",
          paste(bad_cov, collapse = ", "),
          "."
        )
      )
    }
  }
}

validate_engine_packages <- function(engine, packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      paste0(
        "The following packages must be installed to use `engine = '",
        engine,
        "'`: ",
        paste(missing, collapse = ", "),
        "."
      )
    )
  }
}

validate_binary_like <- function(x) {
  if (is.factor(x)) {
    return(nlevels(x) == 2)
  }

  ux <- unique(stats::na.omit(x))
  length(ux) == 2
}

validate_integer_like <- function(x) {
  if (is.factor(x) || is.ordered(x) || is.integer(x)) {
    return(TRUE)
  }

  if (!is.numeric(x)) {
    return(FALSE)
  }

  all(abs(stats::na.omit(x) - round(stats::na.omit(x))) < sqrt(.Machine$double.eps))
}

classify_bart_targets <- function(data, variables) {
  target_data <- data[, variables, drop = FALSE]

  unsupported <- names(Filter(
    function(x) !(is.numeric(x) || is.integer(x) || is.factor(x) || is.ordered(x)),
    target_data
  ))

  if (length(unsupported) > 0) {
    stop(
      paste0(
        "`bart` only supports numeric, integer-like, factor, or ordered synthesized columns. Unsupported columns: ",
        paste(unsupported, collapse = ", "),
        "."
      )
    )
  }

  kinds <- stats::setNames(
    vapply(
      variables,
      function(v) {
        x <- target_data[[v]]

        if (validate_binary_like(x)) {
          return("binary")
        }

        if (is.factor(x) || is.ordered(x)) {
          return("ordinal_approx")
        }

        if (validate_integer_like(x)) {
          return("ordinal_approx")
        }

        "continuous"
      },
      character(1)
    ),
    variables
  )

  numeric_ordinal <- variables[
    kinds == "ordinal_approx" &
      !vapply(target_data, function(x) is.factor(x) || is.ordered(x), logical(1))
  ]

  list(
    kinds = kinds,
    ordinal_approx = variables[kinds == "ordinal_approx"],
    numeric_ordinal = numeric_ordinal
  )
}

validate_factor_data <- function(data, variables) {
  bad <- names(Filter(
    function(x) !(is.factor(x) || is.ordered(x) || validate_integer_like(x)),
    data[, variables, drop = FALSE]
  ))

  if (length(bad) > 0) {
    stop(
      paste0(
        "`factor` currently supports integer-like, factor, or ordered psychometric indicators only. Unsupported columns: ",
        paste(bad, collapse = ", "),
        "."
      )
    )
  }
}

validate_graphical_data <- function(data, variables) {
  bad <- names(Filter(
    function(x) !(is.numeric(x) || is.integer(x) || is.factor(x) || is.ordered(x)),
    data[, variables, drop = FALSE]
  ))

  if (length(bad) > 0) {
    stop(
      paste0(
        "`graphical` currently supports numeric, integer-like, factor, or ordered inputs only. Unsupported columns: ",
        paste(bad, collapse = ", "),
        "."
      )
    )
  }
}

select_predictors <- function(
    target,
    data,
    variables,
    predictor_method,
    predictor_k
) {
  others <- setdiff(variables, target)

  if (length(others) == 0) {
    return(character())
  }

  if (predictor_method == "full") {
    return(others)
  }

  if (predictor_method == "random") {
    return(sample(others, predictor_k, replace = FALSE))
  }

  nmi_values <- normalized_mi_scores(
    target = target,
    data = data,
    variables = variables
  )
  nmi_sorted <- sort(nmi_values, decreasing = TRUE)

  if (predictor_method == "information") {
    return(names(nmi_sorted)[seq_len(min(predictor_k, length(nmi_sorted)))])
  }

  nmi_top <- names(nmi_sorted)[seq_len(min(predictor_k * 5, length(nmi_sorted)))]
  sample(nmi_top, size = min(predictor_k, length(nmi_top)), replace = FALSE)
}

coerce_structural_input <- function(x) {
  if (is.factor(x) || is.ordered(x)) {
    return(as.integer(x))
  }

  if (is.logical(x)) {
    return(as.integer(x))
  }

  as.numeric(x)
}

make_structural_matrix <- function(data, variables) {
  out <- as.data.frame(lapply(data[, variables, drop = FALSE], coerce_structural_input))
  as.matrix(out)
}

make_group_ids <- function(data, covariates = NULL) {
  if (!is.null(covariates)) {
    interaction(data[, covariates], drop = TRUE, lex.order = TRUE)
  } else {
    factor("all")
  }
}

build_group_plan <- function(data, variables, covariates, n_aug, min_group_n = NULL) {
  data$group_id <- make_group_ids(data, covariates)

  if (is.null(min_group_n)) {
    min_group_n <- max(10L, 2L * length(variables))
  }

  proportions <- data |>
    dplyr::group_by(group_id) |>
    dplyr::count(name = "n") |>
    dplyr::ungroup() |>
    dplyr::mutate(
      prop = n / nrow(data),
      target_n = pmax(
        n,
        ceiling(prop * n_aug),
        min_group_n
      ),
      syn = pmax(target_n - n, 0)
    )

  sparse_groups <- proportions$group_id[proportions$n < min_group_n]
  if (length(sparse_groups) > 0) {
    warning(
      paste0(
        "Some covariate-defined subgroups are smaller than the recommended minimum (",
        min_group_n,
        "): ",
        paste(as.character(sparse_groups), collapse = ", "),
        "."
      )
    )
  }

  list(
    data = data,
    groups = unique(as.character(data$group_id)),
    proportions = proportions,
    min_group_n = min_group_n
  )
}

sample_group_feeder <- function(
    group_plan,
    variables,
    noise_method,
    c,
    p,
    f = 0.5,
    r = 1,
    round_int = TRUE,
    noise_bounds = NULL,
    seed = NA
) {
  sampled_data <- purrr::map2_dfr(
    group_plan$groups,
    seq_along(group_plan$groups),
    function(g, i) {
      data_g <- group_plan$data[group_plan$data$group_id == g, , drop = FALSE]
      n_syn_g <- group_plan$proportions$syn[group_plan$proportions$group_id == g]

      if (n_syn_g == 0) {
        return(data_g[0, , drop = FALSE])
      }

      with_seed(seed + i - 1, dplyr::slice_sample(
        .data = data_g,
        n = n_syn_g,
        replace = TRUE
      ))
    }
  )

  noisy_data <- if (nrow(sampled_data) > 0) {
    add_noise(
      data = sampled_data,
      variables = dplyr::all_of(variables),
      method = noise_method,
      c = c,
      p = p,
      f = f,
      r = r,
      round_int = round_int,
      bounds = noise_bounds
    )
  } else {
    sampled_data
  }

  if (nrow(noisy_data) > 0) {
    rbind(group_plan$data, noisy_data)
  } else {
    group_plan$data
  }
}

compute_group_output_sizes <- function(grouped_data, n_aug) {
  group_props <- prop.table(table(grouped_data$group_id))
  group_n <- round(group_props * n_aug)
  diff_n <- n_aug - sum(group_n)
  if (diff_n != 0) {
    adjust_group <- names(which.max(group_n))
    group_n[adjust_group] <- group_n[adjust_group] + diff_n
  }
  group_n
}

make_lightweight_engine_meta <- function(...) {
  list(...)
}

make_predictor_frame <- function(data, columns) {
  if (length(columns) == 0) {
    out <- data.frame(`(Intercept)` = rep(1, nrow(data)))
    return(out)
  }

  data[, columns, drop = FALSE]
}

make_numeric_design_matrix <- function(data) {
  mm <- stats::model.matrix(~ . - 1, data = data)
  if (is.null(dim(mm))) {
    mm <- matrix(mm, ncol = 1)
  }
  mm
}

simulate_latent_gaussian <- function(n, sigma) {
  sigma <- stabilize_correlation_matrix(sigma)

  if (n <= 0) {
    return(matrix(numeric(0), nrow = 0, ncol = ncol(sigma)))
  }

  z <- matrix(stats::rnorm(n * ncol(sigma)), nrow = n, ncol = ncol(sigma))
  z %*% chol(sigma)
}

sample_posterior_column <- function(pp, template, family = NULL) {
  nrow_pp <- nrow(pp)
  ncol_pp <- ncol(pp)
  rows <- sample.int(nrow_pp, ncol_pp, replace = TRUE)
  vec <- pp[cbind(rows, seq_len(ncol_pp))]

  reconstruct_column(vec, template, family = family)
}

reconstruct_column <- function(values, template, family = NULL) {
  observed <- stats::na.omit(template)

  if (is.factor(template)) {
    level_index <- round(values)
    level_index <- pmin(pmax(level_index, 1L), nlevels(template))
    return(factor(
      levels(template)[level_index],
      levels = levels(template),
      ordered = is.ordered(template)
    ))
  }

  if (is.integer(template)) {
    lo <- min(observed, na.rm = TRUE)
    hi <- max(observed, na.rm = TRUE)
    values <- pmin(pmax(values, lo), hi)
    return(as.integer(round(values)))
  }

  is_integer_like <- is.numeric(template) &&
    all(abs(stats::na.omit(template) - round(stats::na.omit(template))) < sqrt(.Machine$double.eps))

  if (is_integer_like) {
    lo <- min(observed, na.rm = TRUE)
    hi <- max(observed, na.rm = TRUE)
    values <- pmin(pmax(values, lo), hi)
    return(as.integer(round(values)))
  }

  if (is.numeric(template) && length(observed) > 0) {
    lo <- min(observed, na.rm = TRUE)
    hi <- max(observed, na.rm = TRUE)
    values <- pmin(pmax(values, lo), hi)
  }

  as.numeric(values)
}

empirical_quantile_map <- function(u, x) {
  observed <- stats::na.omit(x)

  if (length(observed) == 0) {
    stop("Cannot reconstruct empirical margins from an entirely missing column.")
  }

  if (is.factor(x)) {
    sorted_obs <- sort(observed)
    idx <- pmax(1L, ceiling(pmin(pmax(u, 0), 1) * length(sorted_obs)))
    idx[idx == 0] <- 1L
    mapped <- sorted_obs[idx]
    return(factor(
      as.character(mapped),
      levels = levels(x),
      ordered = is.ordered(x)
    ))
  }

  if (is.integer(x)) {
    q <- stats::quantile(
      observed,
      probs = pmin(pmax(u, 0), 1),
      type = 1,
      names = FALSE,
      na.rm = TRUE
    )
    return(as.integer(round(q)))
  }

  is_integer_like <- is.numeric(x) &&
    all(abs(stats::na.omit(x) - round(stats::na.omit(x))) < sqrt(.Machine$double.eps))

  if (is_integer_like) {
    q <- stats::quantile(
      observed,
      probs = pmin(pmax(u, 0), 1),
      type = 1,
      names = FALSE,
      na.rm = TRUE
    )
    return(as.integer(round(q)))
  }

  q <- stats::quantile(
    observed,
    probs = pmin(pmax(u, 0), 1),
    type = 8,
    names = FALSE,
    na.rm = TRUE
  )

  as.numeric(q)
}

stabilize_correlation_matrix <- function(sigma, tol = 1e-8) {
  sigma <- as.matrix(sigma)

  if (length(dim(sigma)) != 2 || nrow(sigma) != ncol(sigma)) {
    stop("`sigma` must be a square matrix.")
  }

  sigma[is.na(sigma)] <- 0
  sigma <- (sigma + t(sigma)) / 2
  diag(sigma) <- 1

  eig <- eigen(sigma, symmetric = TRUE)
  eig$values[eig$values < tol] <- tol

  sigma_pd <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
  sigma_pd <- (sigma_pd + t(sigma_pd)) / 2

  scale_vec <- sqrt(diag(sigma_pd))
  scale_vec[scale_vec < tol] <- 1
  sigma_pd <- diag(1 / scale_vec, nrow = length(scale_vec)) %*%
    sigma_pd %*%
    diag(1 / scale_vec, nrow = length(scale_vec))

  diag(sigma_pd) <- 1

  if (!is.null(colnames(sigma))) {
    colnames(sigma_pd) <- colnames(sigma)
    rownames(sigma_pd) <- rownames(sigma)
  }

  sigma_pd
}

simulate_gaussian_copula <- function(n, sigma) {
  sigma <- stabilize_correlation_matrix(sigma)

  if (n <= 0) {
    return(matrix(numeric(0), nrow = 0, ncol = ncol(sigma)))
  }

  z <- matrix(stats::rnorm(n * ncol(sigma)), nrow = n, ncol = ncol(sigma))
  chol_sigma <- chol(sigma)
  stats::pnorm(z %*% chol_sigma)
}
