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

capture_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    force(expr),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(
    value = value,
    warnings = unique(warnings)
  )
}

quiet_huge_npn <- function(x) {
  out <- NULL
  utils::capture.output(
    out <- suppressMessages(huge::huge.npn(x))
  )
  out
}

is_cscore_composite_list <- function(x) {
  is.list(x) &&
    !is.null(x$lower) &&
    is.list(x$lower) &&
    (is.null(x$higher) || is.list(x$higher))
}

merge_higher_order_blocks <- function(composite_list) {
  lower <- composite_list$lower
  higher <- composite_list$higher

  if (is.null(higher) || length(higher) == 0) {
    return(lower)
  }

  referenced_lower <- unique(unlist(higher, use.names = FALSE))
  missing_lower <- setdiff(referenced_lower, names(lower))
  if (length(missing_lower) > 0) {
    stop(
      paste0(
        "`measurement_blocks$higher` references lower-order composites not present in `measurement_blocks$lower`: ",
        paste(missing_lower, collapse = ", "),
        "."
      )
    )
  }

  higher_blocks <- lapply(higher, function(lower_names) {
    unique(unlist(lower[lower_names], use.names = FALSE))
  })

  remaining_lower <- lower[setdiff(names(lower), referenced_lower)]
  c(remaining_lower, higher_blocks)
}

as_measurement_blocks <- function(measurement_blocks, variables) {
  if (is.null(measurement_blocks)) {
    return(NULL)
  }

  blocks <- if (is_cscore_composite_list(measurement_blocks)) {
    merge_higher_order_blocks(measurement_blocks)
  } else {
    measurement_blocks
  }

  if (!is.list(blocks) || length(blocks) == 0) {
    stop("`measurement_blocks` must be NULL or a non-empty named list.")
  }

  if (is.null(names(blocks)) || any(names(blocks) == "")) {
    stop("`measurement_blocks` must be a named list.")
  }

  blocks <- lapply(blocks, function(x) {
    if (!is.character(x)) {
      stop("Each `measurement_blocks` element must be a character vector.")
    }
    unique(x)
  })

  block_items <- unlist(blocks, use.names = FALSE)
  missing_items <- setdiff(block_items, variables)
  if (length(missing_items) > 0) {
    stop(
      paste0(
        "`measurement_blocks` contains items not present in `variables`: ",
        paste(missing_items, collapse = ", "),
        "."
      )
    )
  }

  duplicated_items <- unique(block_items[duplicated(block_items)])
  if (length(duplicated_items) > 0) {
    stop(
      paste0(
        "`measurement_blocks` items cannot appear in more than one final block: ",
        paste(duplicated_items, collapse = ", "),
        "."
      )
    )
  }

  unblocked <- setdiff(variables, block_items)
  singleton_blocks <- stats::setNames(as.list(unblocked), unblocked)
  c(blocks, singleton_blocks)
}

make_block_shrinkage_matrix <- function(
    variables,
    measurement_blocks,
    between_block_shrink = 0.5
) {
  if (
    !is.numeric(between_block_shrink) ||
      length(between_block_shrink) != 1 ||
      is.na(between_block_shrink) ||
      between_block_shrink < 0 ||
      between_block_shrink > 1
  ) {
    stop("`between_block_shrink` must be a single numeric value in [0, 1].")
  }

  block_lookup <- stats::setNames(rep(NA_character_, length(variables)), variables)
  for (block_name in names(measurement_blocks)) {
    block_lookup[measurement_blocks[[block_name]]] <- block_name
  }

  shrinkage <- matrix(
    between_block_shrink,
    nrow = length(variables),
    ncol = length(variables),
    dimnames = list(variables, variables)
  )
  same_block <- outer(block_lookup, block_lookup, "==")
  shrinkage[same_block] <- 1
  diag(shrinkage) <- 1
  shrinkage
}

apply_block_shrinkage <- function(cor_mat, shrinkage_matrix) {
  shrinkage_matrix <- shrinkage_matrix[rownames(cor_mat), colnames(cor_mat), drop = FALSE]
  out <- cor_mat * shrinkage_matrix
  diag(out) <- 1
  out <- (out + t(out)) / 2
  stabilize_correlation_matrix(out)
}

restore_within_block_correlations <- function(cor_mat, reference_cor, measurement_blocks) {
  if (is.null(measurement_blocks)) {
    return(cor_mat)
  }

  out <- as.matrix(cor_mat)
  reference_cor <- as.matrix(reference_cor)
  reference_cor <- reference_cor[rownames(out), colnames(out), drop = FALSE]

  for (block in measurement_blocks) {
    block <- intersect(block, rownames(out))
    if (length(block) >= 2) {
      out[block, block] <- reference_cor[block, block, drop = FALSE]
    }
  }

  diag(out) <- 1
  out <- (out + t(out)) / 2
  out
}

estimate_factor_item_correlation <- function(data, variables, correlation_method = "polychoric") {
  struct_active <- as.data.frame(lapply(data[, variables, drop = FALSE], coerce_structural_input))

  cor_mat <- if (identical(correlation_method, "polychoric")) {
    tryCatch(
      psych::polychoric(struct_active)$rho,
      error = function(e) stats::cor(struct_active, use = "pairwise.complete.obs")
    )
  } else {
    stats::cor(struct_active, use = "pairwise.complete.obs")
  }

  cor_mat <- as.matrix(cor_mat)
  rownames(cor_mat) <- variables
  colnames(cor_mat) <- variables
  stabilize_correlation_matrix(cor_mat)
}

make_measurement_block_scores <- function(data, measurement_blocks) {
  scores <- lapply(measurement_blocks, function(block) {
    struct_block <- as.data.frame(lapply(data[, block, drop = FALSE], coerce_structural_input))
    if (length(block) == 1) {
      raw_score <- struct_block[[1]]
    } else {
      z_block <- as.data.frame(lapply(struct_block, function(x) {
        x_sd <- stats::sd(x, na.rm = TRUE)
        if (!is.finite(x_sd) || x_sd <= sqrt(.Machine$double.eps)) {
          return(rep(NA_real_, length(x)))
        }
        as.numeric(scale(x))
      }))
      raw_score <- rowMeans(z_block, na.rm = TRUE)
      raw_score[!is.finite(raw_score)] <- NA_real_
    }

    score_sd <- stats::sd(raw_score, na.rm = TRUE)
    if (!is.finite(score_sd) || score_sd <= sqrt(.Machine$double.eps)) {
      return(rep(0, length(raw_score)))
    }
    as.numeric(scale(raw_score))
  })

  out <- as.data.frame(scores)
  names(out) <- names(measurement_blocks)
  out
}

apply_between_block_shrinkage <- function(cor_mat, between_block_shrink = 0.5) {
  if (
    !is.numeric(between_block_shrink) ||
      length(between_block_shrink) != 1 ||
      is.na(between_block_shrink) ||
      between_block_shrink < 0 ||
      between_block_shrink > 1
  ) {
    stop("`between_block_shrink` must be a single numeric value in [0, 1].")
  }

  out <- as.matrix(cor_mat)
  out[is.na(out)] <- 0
  out <- (out + t(out)) / 2
  diag(out) <- 1
  out[row(out) != col(out)] <- out[row(out) != col(out)] * between_block_shrink
  stabilize_correlation_matrix(out)
}

make_brm_measurement_plan <- function(measurement_blocks, variables) {
  if (is.null(measurement_blocks)) {
    return(NULL)
  }

  if (is_cscore_composite_list(measurement_blocks)) {
    return(make_brm_cscore_measurement_plan(measurement_blocks, variables))
  }

  if (!is.list(measurement_blocks) || length(measurement_blocks) == 0) {
    stop("`measurement_blocks` must be NULL or a non-empty named list.")
  }
  if (is.null(names(measurement_blocks)) || any(names(measurement_blocks) == "")) {
    stop("`measurement_blocks` must be a named list.")
  }

  leaf_blocks <- lapply(measurement_blocks, function(items) {
    if (!is.character(items)) {
      stop("Each `measurement_blocks` element must be a character vector.")
    }
    unique(items)
  })
  validate_brm_leaf_blocks(leaf_blocks, variables)
  leaf_blocks <- add_brm_singleton_leaf_blocks(leaf_blocks, variables)

  structural_families <- stats::setNames(
    lapply(names(leaf_blocks), function(block_name) {
      list(
        leaf_blocks = block_name,
        variables = leaf_blocks[[block_name]],
        source = "plain"
      )
    }),
    names(leaf_blocks)
  )

  make_brm_measurement_plan_object(
    leaf_blocks = leaf_blocks,
    structural_families = structural_families,
    ignored_composites = character(),
    input_type = "plain"
  )
}

make_measurement_score_plan <- function(measurement_blocks, variables) {
  plan <- make_brm_measurement_plan(
    measurement_blocks = measurement_blocks,
    variables = variables
  )

  if (!is.null(plan) && identical(plan$input_type, "cscore")) {
    sources <- vapply(
      plan$structural_families,
      function(x) x$source,
      character(1)
    )
    family_order <- c(
      names(sources)[sources == "lower"],
      names(sources)[sources == "higher"],
      names(sources)[sources == "singleton"]
    )
    plan$structural_families <- plan$structural_families[family_order]
    plan$family_names <- family_order
  }

  plan
}

make_brm_cscore_measurement_plan <- function(composite_list, variables) {
  lower <- composite_list$lower
  higher <- composite_list$higher

  if (!is.list(lower) || length(lower) == 0) {
    stop("`measurement_blocks$lower` must be a non-empty list.")
  }
  if (is.null(names(lower)) || any(names(lower) == "")) {
    stop("`measurement_blocks$lower` must be a named list.")
  }

  ignored <- character()
  leaf_blocks <- list()
  for (block_name in names(lower)) {
    items <- lower[[block_name]]
    if (!is.character(items)) {
      stop("Each `measurement_blocks$lower` element must be a character vector.")
    }
    items <- unique(items)
    in_variables <- items %in% variables

    if (all(in_variables)) {
      leaf_blocks[[block_name]] <- items
    } else if (!any(in_variables)) {
      ignored <- c(ignored, block_name)
    } else {
      stop(
        paste0(
          "`measurement_blocks$lower$", block_name,
          "` mixes synthesized and non-synthesized items. Either include all items in `variables` or remove the composite."
        )
      )
    }
  }

  validate_brm_leaf_blocks(leaf_blocks, variables)

  structural_families <- list()
  referenced_valid_lower <- character()
  if (!is.null(higher) && length(higher) > 0) {
    if (!is.list(higher) || is.null(names(higher)) || any(names(higher) == "")) {
      stop("`measurement_blocks$higher` must be a named list when supplied.")
    }

    for (family_name in names(higher)) {
      lower_refs <- higher[[family_name]]
      if (!is.character(lower_refs)) {
        stop("Each `measurement_blocks$higher` element must be a character vector.")
      }
      lower_refs <- unique(lower_refs)

      missing_lower <- setdiff(lower_refs, names(lower))
      if (length(missing_lower) > 0) {
        stop(
          paste0(
            "`measurement_blocks$higher` references lower-order composites not present in `measurement_blocks$lower`: ",
            paste(missing_lower, collapse = ", "),
            "."
          )
        )
      }

      valid_refs <- intersect(lower_refs, names(leaf_blocks))
      if (length(valid_refs) == 0) {
        ignored <- c(ignored, family_name)
        next
      }

      structural_families[[family_name]] <- list(
        leaf_blocks = valid_refs,
        variables = unique(unlist(leaf_blocks[valid_refs], use.names = FALSE)),
        source = "higher"
      )
      referenced_valid_lower <- union(referenced_valid_lower, valid_refs)
    }
  }

  remaining_lower <- setdiff(names(leaf_blocks), referenced_valid_lower)
  for (block_name in remaining_lower) {
    structural_families[[block_name]] <- list(
      leaf_blocks = block_name,
      variables = leaf_blocks[[block_name]],
      source = "lower"
    )
  }

  leaf_blocks <- add_brm_singleton_leaf_blocks(leaf_blocks, variables)
  singleton_names <- setdiff(names(leaf_blocks), c(names(lower), names(structural_families)))
  for (block_name in singleton_names) {
    structural_families[[block_name]] <- list(
      leaf_blocks = block_name,
      variables = leaf_blocks[[block_name]],
      source = "singleton"
    )
  }

  make_brm_measurement_plan_object(
    leaf_blocks = leaf_blocks,
    structural_families = structural_families,
    ignored_composites = unique(ignored),
    input_type = "cscore"
  )
}

validate_brm_leaf_blocks <- function(leaf_blocks, variables) {
  if (length(leaf_blocks) == 0) {
    return(invisible(TRUE))
  }

  block_items <- unlist(leaf_blocks, use.names = FALSE)
  missing_items <- setdiff(block_items, variables)
  if (length(missing_items) > 0) {
    stop(
      paste0(
        "`measurement_blocks` contains items not present in `variables`: ",
        paste(missing_items, collapse = ", "),
        "."
      )
    )
  }

  duplicated_items <- unique(block_items[duplicated(block_items)])
  if (length(duplicated_items) > 0) {
    stop(
      paste0(
        "`measurement_blocks` items cannot appear in more than one measurement leaf block: ",
        paste(duplicated_items, collapse = ", "),
        "."
      )
    )
  }

  invisible(TRUE)
}

add_brm_singleton_leaf_blocks <- function(leaf_blocks, variables) {
  block_items <- unlist(leaf_blocks, use.names = FALSE)
  unblocked <- setdiff(variables, block_items)
  if (length(unblocked) == 0) {
    return(leaf_blocks)
  }

  c(leaf_blocks, stats::setNames(as.list(unblocked), unblocked))
}

make_brm_measurement_plan_object <- function(
    leaf_blocks,
    structural_families,
    ignored_composites,
    input_type
) {
  if (length(leaf_blocks) == 0) {
    stop("No synthesized variables were found in `measurement_blocks`.")
  }
  if (length(structural_families) == 0) {
    stop("No structural families could be formed from `measurement_blocks`.")
  }

  family_leaves <- unlist(
    lapply(structural_families, function(x) x$leaf_blocks),
    use.names = FALSE
  )
  duplicated_family_leaves <- unique(family_leaves[duplicated(family_leaves)])
  if (length(duplicated_family_leaves) > 0) {
    stop(
      paste0(
        "Each lower-order measurement block can belong to only one final structural family: ",
        paste(duplicated_family_leaves, collapse = ", "),
        "."
      )
    )
  }

  leaf_to_family <- stats::setNames(rep(NA_character_, length(leaf_blocks)), names(leaf_blocks))
  for (family_name in names(structural_families)) {
    leaf_to_family[structural_families[[family_name]]$leaf_blocks] <- family_name
  }

  list(
    leaf_blocks = leaf_blocks,
    structural_families = structural_families,
    family_names = names(structural_families),
    leaf_to_family = leaf_to_family,
    ignored_composites = ignored_composites,
    input_type = input_type
  )
}

make_brm_leaf_score_frame <- function(data, leaf_blocks) {
  make_measurement_block_scores(data = data, measurement_blocks = leaf_blocks)
}

make_brm_structural_score_frame <- function(leaf_scores, measurement_plan) {
  family_scores <- lapply(measurement_plan$structural_families, function(family) {
    leaves <- family$leaf_blocks
    if (length(leaves) == 1) {
      raw_score <- leaf_scores[[leaves]]
    } else {
      raw_score <- rowMeans(leaf_scores[, leaves, drop = FALSE], na.rm = TRUE)
      raw_score[!is.finite(raw_score)] <- NA_real_
    }

    standardize_brm_score(raw_score)
  })

  out <- as.data.frame(family_scores)
  names(out) <- names(measurement_plan$structural_families)
  out
}

make_measurement_leaf_score_frame <- function(data, measurement_plan) {
  make_brm_leaf_score_frame(
    data = data,
    leaf_blocks = measurement_plan$leaf_blocks
  )
}

make_measurement_structural_score_frame <- function(data, measurement_plan) {
  leaf_scores <- make_measurement_leaf_score_frame(
    data = data,
    measurement_plan = measurement_plan
  )
  make_brm_structural_score_frame(
    leaf_scores = leaf_scores,
    measurement_plan = measurement_plan
  )
}

standardize_brm_score <- function(x) {
  x <- as.numeric(x)
  x_mean <- mean(x, na.rm = TRUE)
  x_sd <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(x_sd) || x_sd <= sqrt(.Machine$double.eps)) {
    return(rep(0, length(x)))
  }
  out <- (x - x_mean) / x_sd
  out[!is.finite(out)] <- 0
  out
}

make_brm_score_column_map <- function(family_names) {
  safe_names <- make.names(family_names, unique = TRUE)
  stats::setNames(paste0("sydni_brm_score_", safe_names), family_names)
}

select_brm_score_predictors <- function(
    target,
    score_data,
    score_columns,
    predictor_method,
    predictor_k,
    between_block_shrink
) {
  candidates <- setdiff(score_columns, target)
  between_k <- min(
    length(candidates),
    max(0L, as.integer(floor(predictor_k * between_block_shrink)))
  )

  predictors <- if (between_k == 0 || length(candidates) == 0) {
    character()
  } else if (identical(predictor_method, "full")) {
    candidates[seq_len(between_k)]
  } else {
    select_predictors(
      target = target,
      data = score_data,
      variables = c(target, candidates),
      predictor_method = predictor_method,
      predictor_k = between_k
    )
  }

  list(
    predictors = predictors,
    within_block_predictors = character(),
    between_block_predictors = predictors,
    between_k = between_k,
    target_block = target
  )
}

make_brm_structural_residual_donors <- function(
    n_obs,
    n_aug,
    family_names,
    seed = NA
) {
  if (n_obs <= 0) {
    stop("`n_obs` must be positive.")
  }
  if (n_aug < 0) {
    stop("`n_aug` must be non-negative.")
  }

  stats::setNames(
    lapply(seq_along(family_names), function(i) {
      with_seed(seed + i - 1L, {
        sample.int(n_obs, size = n_aug, replace = TRUE)
      })
    }),
    family_names
  )
}

estimate_brm_score_correlation <- function(score_data) {
  score_data <- as.data.frame(score_data)
  p <- ncol(score_data)
  if (p == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  if (p == 1) {
    out <- matrix(1, nrow = 1, ncol = 1)
    rownames(out) <- colnames(out) <- names(score_data)
    return(out)
  }

  cor_mat <- suppressWarnings(stats::cor(
    score_data,
    use = "pairwise.complete.obs"
  ))
  cor_mat <- as.matrix(cor_mat)
  cor_mat[!is.finite(cor_mat)] <- 0
  cor_mat <- (cor_mat + t(cor_mat)) / 2
  diag(cor_mat) <- 1
  rownames(cor_mat) <- colnames(cor_mat) <- names(score_data)
  stabilize_correlation_matrix(cor_mat)
}

calibrate_brm_structural_score_correlations <- function(
    score_data,
    target_cor,
    seed = NA,
    calibration_draws = 25L
) {
  score_data <- as.data.frame(score_data)
  p <- ncol(score_data)
  n <- nrow(score_data)
  calibration_draws <- max(1L, as.integer(calibration_draws))
  target_cor <- target_cor[names(score_data), names(score_data), drop = FALSE]
  target_cor <- stabilize_correlation_matrix(target_cor)
  raw_cor <- estimate_brm_score_correlation(score_data)
  raw_error <- summarize_brm_structural_correlation_error(raw_cor, target_cor)

  if (p < 2 || n < 3) {
    attr(score_data, "calibration") <- list(
      correlation = raw_cor,
      error = raw_error$mean,
      max_error = raw_error$max,
      raw_error = raw_error$mean,
      raw_max_error = raw_error$max,
      draws = calibration_draws,
      selected_draw = NA_integer_
    )
    return(score_data)
  }

  best <- NULL

  for (draw in seq_len(calibration_draws)) {
    draw_seed <- if (length(seed) == 1 && is.finite(seed) && !is.na(seed)) {
      seed + draw - 1L
    } else {
      NA_real_
    }

    candidate <- calibrate_brm_structural_score_correlations_once(
      score_data = score_data,
      target_cor = target_cor,
      seed = draw_seed
    )
    candidate_cor <- estimate_brm_score_correlation(candidate)
    candidate_error <- summarize_brm_structural_correlation_error(candidate_cor, target_cor)

    if (
      is.null(best) ||
        candidate_error$max < best$max_error ||
        (
          isTRUE(all.equal(candidate_error$max, best$max_error, tolerance = 1e-12)) &&
            candidate_error$mean < best$error
        )
    ) {
      best <- list(
        data = candidate,
        correlation = candidate_cor,
        error = candidate_error$mean,
        max_error = candidate_error$max,
        selected_draw = draw
      )
    }
  }

  attr(best$data, "calibration") <- list(
    correlation = best$correlation,
    error = best$error,
    max_error = best$max_error,
    raw_error = raw_error$mean,
    raw_max_error = raw_error$max,
    draws = calibration_draws,
    selected_draw = best$selected_draw
  )
  best$data
}

calibrate_brm_structural_score_correlations_once <- function(
    score_data,
    target_cor,
    seed = NA
) {
  n <- nrow(score_data)
  p <- ncol(score_data)

  latent <- with_seed(seed, {
    matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  })
  latent <- scale(latent)
  latent[!is.finite(latent)] <- 0

  latent_cor <- suppressWarnings(stats::cor(latent))
  latent_cor[!is.finite(latent_cor)] <- 0
  latent_cor <- stabilize_correlation_matrix(latent_cor)
  latent_chol <- chol(latent_cor)
  target_chol <- chol(target_cor)
  latent_target <- latent %*% solve(latent_chol, target_chol)
  colnames(latent_target) <- names(score_data)

  out <- score_data
  for (col in names(score_data)) {
    target_order <- order(latent_target[, col], seq_len(n), na.last = TRUE)
    values <- sort(score_data[[col]], na.last = TRUE)
    calibrated <- score_data[[col]]
    calibrated[target_order] <- values
    out[[col]] <- calibrated
  }

  out
}

summarize_brm_structural_correlation_error <- function(cor_mat, target_cor) {
  cor_mat <- cor_mat[rownames(target_cor), colnames(target_cor), drop = FALSE]
  err <- abs(cor_mat - target_cor)
  err[!is.finite(err)] <- 0
  diag(err) <- 0
  offdiag <- err[row(err) != col(err)]

  if (length(offdiag) == 0) {
    return(list(mean = 0, max = 0))
  }

  list(
    mean = mean(offdiag),
    max = max(offdiag)
  )
}

make_measurement_score_permutation <- function(source_scores, target_scores) {
  if (length(source_scores) != length(target_scores)) {
    stop("`source_scores` and `target_scores` must have the same length.")
  }

  n <- length(source_scores)
  if (n == 0) {
    return(integer())
  }

  source_scores <- as.numeric(source_scores)
  target_scores <- as.numeric(target_scores)
  source_scores[!is.finite(source_scores)] <- 0
  target_scores[!is.finite(target_scores)] <- 0

  source_order <- order(source_scores, seq_len(n), na.last = TRUE)
  target_order <- order(target_scores, seq_len(n), na.last = TRUE)
  permutation <- integer(n)
  permutation[target_order] <- source_order
  permutation
}

apply_measurement_structural_bundle_permutation <- function(
    data,
    source_scores,
    target_scores,
    measurement_plan
) {
  out <- data
  permutations <- list()

  for (family_name in measurement_plan$family_names) {
    family_vars <- measurement_plan$structural_families[[family_name]]$variables
    permutation <- make_measurement_score_permutation(
      source_scores = source_scores[[family_name]],
      target_scores = target_scores[[family_name]]
    )
    out[, family_vars] <- data[permutation, family_vars, drop = FALSE]
    permutations[[family_name]] <- permutation
  }

  attr(out, "family_permutations") <- permutations
  out
}

estimate_measurement_within_family_correlation_change <- function(
    before,
    after,
    measurement_plan
) {
  changes <- vapply(measurement_plan$family_names, function(family_name) {
    family_vars <- measurement_plan$structural_families[[family_name]]$variables
    if (length(family_vars) < 2) {
      return(0)
    }

    before_frame <- as.data.frame(lapply(
      before[, family_vars, drop = FALSE],
      coerce_structural_input
    ))
    after_frame <- as.data.frame(lapply(
      after[, family_vars, drop = FALSE],
      coerce_structural_input
    ))
    before_cor <- suppressWarnings(stats::cor(before_frame, use = "pairwise.complete.obs"))
    after_cor <- suppressWarnings(stats::cor(after_frame, use = "pairwise.complete.obs"))
    differences <- abs(before_cor - after_cor)
    differences <- differences[is.finite(differences)]
    if (length(differences) == 0) 0 else max(differences)
  }, numeric(1))

  list(
    by_family = changes,
    max = if (length(changes) == 0) 0 else max(changes)
  )
}

estimate_measurement_score_covariate_correlations <- function(
    score_data,
    data,
    covariates
) {
  score_data <- as.data.frame(score_data)
  if (length(covariates) == 0) {
    return(matrix(
      numeric(0),
      nrow = ncol(score_data),
      ncol = 0,
      dimnames = list(names(score_data), character())
    ))
  }

  covariate_frame <- data[, covariates, drop = FALSE]
  covariate_matrix <- tryCatch(
    make_numeric_design_matrix(covariate_frame),
    error = function(e) NULL
  )
  if (is.null(covariate_matrix) || nrow(covariate_matrix) != nrow(score_data)) {
    return(matrix(
      NA_real_,
      nrow = ncol(score_data),
      ncol = 0,
      dimnames = list(names(score_data), character())
    ))
  }

  colnames(covariate_matrix) <- make.unique(colnames(covariate_matrix))
  combined <- cbind(as.matrix(score_data), covariate_matrix)
  correlations <- suppressWarnings(stats::cor(combined, use = "pairwise.complete.obs"))
  correlations[
    names(score_data),
    colnames(covariate_matrix),
    drop = FALSE
  ]
}

calibrate_measurement_structural_bundles <- function(
    data,
    measurement_plan,
    target_cor,
    seed = NA,
    calibration_draws = 25L
) {
  calibration_draws <- max(1L, as.integer(calibration_draws))
  raw_scores <- make_measurement_structural_score_frame(
    data = data,
    measurement_plan = measurement_plan
  )
  target_cor <- target_cor[
    measurement_plan$family_names,
    measurement_plan$family_names,
    drop = FALSE
  ]
  target_cor <- stabilize_correlation_matrix(target_cor)
  raw_cor <- estimate_brm_score_correlation(raw_scores)
  raw_error <- summarize_brm_structural_correlation_error(raw_cor, target_cor)

  best <- list(
    data = data,
    scores = raw_scores,
    correlation = raw_cor,
    error = raw_error$mean,
    max_error = raw_error$max,
    selected_draw = 0L,
    permutations = stats::setNames(
      rep(list(seq_len(nrow(data))), length(measurement_plan$family_names)),
      measurement_plan$family_names
    )
  )

  fallback_reason <- NULL
  if (ncol(raw_scores) < 2) {
    fallback_reason <- "fewer_than_two_structural_families"
  } else if (nrow(raw_scores) < 3) {
    fallback_reason <- "fewer_than_three_synthetic_rows"
  } else {
    for (draw in seq_len(calibration_draws)) {
      draw_seed <- if (length(seed) == 1 && is.finite(seed) && !is.na(seed)) {
        seed + draw - 1L
      } else {
        NA_real_
      }
      candidate_score_target <- calibrate_brm_structural_score_correlations_once(
        score_data = raw_scores,
        target_cor = target_cor,
        seed = draw_seed
      )
      candidate_data <- apply_measurement_structural_bundle_permutation(
        data = data,
        source_scores = raw_scores,
        target_scores = candidate_score_target,
        measurement_plan = measurement_plan
      )
      permutations <- attr(candidate_data, "family_permutations")
      attr(candidate_data, "family_permutations") <- NULL
      candidate_scores <- make_measurement_structural_score_frame(
        data = candidate_data,
        measurement_plan = measurement_plan
      )
      candidate_cor <- estimate_brm_score_correlation(candidate_scores)
      candidate_error <- summarize_brm_structural_correlation_error(
        candidate_cor,
        target_cor
      )

      if (
        candidate_error$max < best$max_error ||
          (
            isTRUE(all.equal(candidate_error$max, best$max_error, tolerance = 1e-12)) &&
              candidate_error$mean < best$error
          )
      ) {
        best <- list(
          data = candidate_data,
          scores = candidate_scores,
          correlation = candidate_cor,
          error = candidate_error$mean,
          max_error = candidate_error$max,
          selected_draw = draw,
          permutations = permutations
        )
      }
    }

    if (identical(best$selected_draw, 0L)) {
      fallback_reason <- "raw_candidate_selected"
    }
  }

  within_change <- estimate_measurement_within_family_correlation_change(
    before = data,
    after = best$data,
    measurement_plan = measurement_plan
  )
  attr(best$data, "calibration") <- list(
    scores_raw = raw_scores,
    scores = best$scores,
    correlation_raw = raw_cor,
    correlation = best$correlation,
    error_raw = raw_error$mean,
    error = best$error,
    max_error_raw = raw_error$max,
    max_error = best$max_error,
    draws = calibration_draws,
    selected_draw = best$selected_draw,
    fallback_reason = fallback_reason,
    permutations = best$permutations,
    within_family_correlation_change = within_change$by_family,
    within_family_correlation_max_change = within_change$max
  )
  best$data
}

simulate_brm_leaf_scores <- function(
    synthetic_family_scores,
    observed_family_scores,
    observed_leaf_scores,
    measurement_plan,
    seed = NA
) {
  n <- nrow(synthetic_family_scores)
  out <- as.data.frame(matrix(
    NA_real_,
    nrow = n,
    ncol = length(measurement_plan$leaf_blocks)
  ))
  names(out) <- names(measurement_plan$leaf_blocks)

  for (family_index in seq_along(measurement_plan$structural_families)) {
    family_name <- names(measurement_plan$structural_families)[[family_index]]
    leaves <- measurement_plan$structural_families[[family_name]]$leaf_blocks
    family_score <- synthetic_family_scores[[family_name]]

    if (length(leaves) == 1) {
      out[[leaves]] <- family_score
      next
    }

    observed_block <- data.frame(
      family = observed_family_scores[[family_name]],
      observed_leaf_scores[, leaves, drop = FALSE],
      check.names = FALSE
    )
    means <- vapply(observed_block, mean, numeric(1), na.rm = TRUE)
    sigma <- stats::cov(observed_block, use = "pairwise.complete.obs")
    sigma <- as.matrix(sigma)
    sigma[!is.finite(sigma)] <- 0
    sigma <- stabilize_covariance_matrix(sigma)

    sigma_ff <- max(sigma["family", "family"], 1e-8)
    sigma_lf <- sigma[leaves, "family", drop = FALSE]
    sigma_ll <- sigma[leaves, leaves, drop = FALSE]
    cond_mean <- matrix(family_score - means[["family"]], ncol = 1) %*%
      t(sigma_lf / sigma_ff)
    cond_mean <- sweep(cond_mean, 2, means[leaves], "+")
    cond_cov <- sigma_ll - sigma_lf %*% t(sigma_lf) / sigma_ff
    cond_cov <- stabilize_covariance_matrix(cond_cov)

    residuals <- with_seed(seed + family_index, {
      matrix(stats::rnorm(n * length(leaves)), nrow = n, ncol = length(leaves))
    }) %*% chol(cond_cov)

    out[, leaves] <- cond_mean + residuals
  }

  out
}

reconstruct_brm_items_from_leaf_scores <- function(
    synthetic_leaf_scores,
    observed_leaf_scores,
    data,
    leaf_blocks,
    variables,
    seed = NA,
    donor_k = 5L
) {
  out <- data[rep(1, nrow(synthetic_leaf_scores)), variables, drop = FALSE]
  donor_rows <- list()

  for (block_index in seq_along(leaf_blocks)) {
    block_name <- names(leaf_blocks)[[block_index]]
    block_vars <- leaf_blocks[[block_name]]
    target_scores <- synthetic_leaf_scores[[block_name]]
    observed_scores <- observed_leaf_scores[[block_name]]
    valid_observed <- which(is.finite(observed_scores))
    if (length(valid_observed) == 0) {
      valid_observed <- seq_len(nrow(data))
      observed_scores <- rep(0, nrow(data))
    }

    donors <- if (length(block_vars) == 2L) {
      match_two_item_brm_donors_by_rank(
        target_scores = target_scores,
        observed_scores = observed_scores,
        valid_observed = valid_observed,
        seed = seed + block_index
      )
    } else {
      with_seed(seed + block_index, {
        vapply(seq_along(target_scores), function(i) {
          target_value <- target_scores[[i]]
          if (!is.finite(target_value)) {
            target_value <- 0
          }
          distances <- abs(observed_scores[valid_observed] - target_value)
          nearest <- valid_observed[order(distances, na.last = NA)]
          nearest <- head(nearest, max(1L, min(as.integer(donor_k), length(nearest))))
          sample(nearest, size = 1L)
        }, integer(1))
      })
    }

    out[, block_vars] <- data[donors, block_vars, drop = FALSE]
    donor_rows[[block_name]] <- donors
  }

  attr(out, "donor_rows") <- donor_rows
  out
}

match_two_item_brm_donors_by_rank <- function(
    target_scores,
    observed_scores,
    valid_observed,
    seed = NA
) {
  n_target <- length(target_scores)
  if (n_target == 0) {
    return(integer())
  }

  observed_scores_valid <- observed_scores[valid_observed]
  observed_scores_valid[!is.finite(observed_scores_valid)] <- 0
  target_scores_valid <- target_scores
  target_scores_valid[!is.finite(target_scores_valid)] <- 0

  matched <- with_seed(seed, {
    observed_order <- valid_observed[order(
      observed_scores_valid,
      stats::runif(length(valid_observed)),
      na.last = TRUE
    )]
    target_rank <- rank(target_scores_valid, ties.method = "random")
    probs <- (target_rank - 0.5) / n_target
    donor_index <- pmin(
      length(observed_order),
      pmax(1L, ceiling(probs * length(observed_order)))
    )
    observed_order[donor_index]
  })

  matched
}

stabilize_covariance_matrix <- function(sigma, tol = 1e-8) {
  sigma <- as.matrix(sigma)
  sigma[!is.finite(sigma)] <- 0
  sigma <- (sigma + t(sigma)) / 2

  eig <- eigen(sigma, symmetric = TRUE)
  eig$values[eig$values < tol] <- tol
  out <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
  out <- (out + t(out)) / 2
  rownames(out) <- rownames(sigma)
  colnames(out) <- colnames(sigma)
  out
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

select_measurement_predictors <- function(
    target,
    data,
    variables,
    predictor_method,
    predictor_k,
    measurement_blocks,
    between_block_shrink
) {
  target_block <- names(Filter(function(block) target %in% block, measurement_blocks))
  if (length(target_block) == 0) {
    within_block_predictors <- character()
  } else {
    within_block_predictors <- setdiff(measurement_blocks[[target_block[[1]]]], target)
  }

  within_block_predictors <- intersect(within_block_predictors, variables)
  between_candidates <- setdiff(variables, c(target, within_block_predictors))
  between_k <- min(
    length(between_candidates),
    max(0L, as.integer(floor(predictor_k * between_block_shrink)))
  )

  between_block_predictors <- if (between_k == 0 || length(between_candidates) == 0) {
    character()
  } else if (identical(predictor_method, "full")) {
    between_candidates[seq_len(between_k)]
  } else {
    select_predictors(
      target = target,
      data = data,
      variables = c(target, between_candidates),
      predictor_method = predictor_method,
      predictor_k = between_k
    )
  }

  list(
    predictors = unique(c(within_block_predictors, between_block_predictors)),
    within_block_predictors = within_block_predictors,
    between_block_predictors = between_block_predictors,
    between_k = between_k,
    target_block = if (length(target_block) == 0) NA_character_ else target_block[[1]]
  )
}

is_brm_coefficient_prior_safe <- function(data, term) {
  if (!term %in% names(data)) {
    return(FALSE)
  }

  x <- data[[term]]
  (is.numeric(x) || is.integer(x) || is.logical(x)) && identical(make.names(term), term)
}

make_brm_role_priors <- function(
    data,
    predictor_roles,
    covariates,
    base_sd,
    shrink_sd,
    between_block_shrink
) {
  rhs_terms <- c(
    predictor_roles$within_block_predictors,
    predictor_roles$between_block_predictors,
    covariates
  )
  rhs_terms <- unique(rhs_terms)

  unsafe_terms <- rhs_terms[
    !vapply(rhs_terms, function(term) is_brm_coefficient_prior_safe(data, term), logical(1))
  ]

  if (length(unsafe_terms) > 0) {
    return(list(
      priors = do.call(brms::prior, list(
        prior = bquote(normal(0, .(shrink_sd))),
        class = "b"
      )),
      meta = list(
        prior_mode = "global_fallback",
        same_block_sd = NA_real_,
        between_block_sd = NA_real_,
        covariate_sd = NA_real_,
        fallback_terms = unsafe_terms
      )
    ))
  }

  same_block_sd <- base_sd
  between_block_sd <- shrink_sd * between_block_shrink
  covariate_sd <- shrink_sd

  priors <- list()
  for (term in predictor_roles$within_block_predictors) {
    priors[[length(priors) + 1L]] <- do.call(brms::prior, list(
      prior = bquote(normal(0, .(same_block_sd))),
      class = "b",
      coef = term
    ))
  }
  for (term in predictor_roles$between_block_predictors) {
    priors[[length(priors) + 1L]] <- do.call(brms::prior, list(
      prior = bquote(normal(0, .(between_block_sd))),
      class = "b",
      coef = term
    ))
  }
  for (term in covariates) {
    priors[[length(priors) + 1L]] <- do.call(brms::prior, list(
      prior = bquote(normal(0, .(covariate_sd))),
      class = "b",
      coef = term
    ))
  }

  list(
    priors = if (length(priors) == 0) NULL else do.call(c, priors),
    meta = list(
      prior_mode = "role_specific",
      same_block_sd = same_block_sd,
      between_block_sd = between_block_sd,
      covariate_sd = covariate_sd,
      fallback_terms = character()
    )
  )
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
    sparse_groups_chr <- as.character(sparse_groups)
    shown_groups <- head(sparse_groups_chr, 10)
    hidden_n <- length(sparse_groups_chr) - length(shown_groups)
    hidden_text <- if (hidden_n > 0) {
      paste0(", and ", hidden_n, " more")
    } else {
      ""
    }

    warning(
      paste0(
        "Some covariate-defined subgroups are smaller than the recommended minimum (",
        min_group_n,
        "). Showing ",
        length(shown_groups),
        " of ",
        length(sparse_groups_chr),
        " sparse subgroup(s): ",
        paste(shown_groups, collapse = ", "),
        hidden_text,
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

sample_coherent_gaussian_epred_column <- function(
    epred_new,
    epred_train,
    y_train,
    template,
    donor_rows,
    epred_rows = NULL,
    residual_scale = 1
) {
  if (length(dim(epred_new)) != 2 || length(dim(epred_train)) != 2) {
    stop("`epred_new` and `epred_train` must be posterior expectation matrices.")
  }

  rows <- epred_rows
  if (is.null(rows)) {
    rows <- sample.int(nrow(epred_new), ncol(epred_new), replace = TRUE)
  }
  if (length(rows) != ncol(epred_new)) {
    stop("`epred_rows` must have one entry per posterior expectation column.")
  }
  if (any(is.na(rows)) || any(rows < 1) || any(rows > nrow(epred_new))) {
    stop("`epred_rows` contains invalid posterior expectation row indices.")
  }
  raw <- epred_new[cbind(rows, seq_len(ncol(epred_new)))]

  fitted_mean <- colMeans(epred_train, na.rm = TRUE)
  residuals <- as.numeric(y_train) - fitted_mean
  residuals[!is.finite(residuals)] <- 0

  if (length(residuals) > 0 && stats::sd(residuals) > sqrt(.Machine$double.eps)) {
    residuals <- residuals - mean(residuals)
    raw <- raw + residuals[donor_rows] * residual_scale
  }

  reconstruct_column(raw, template)
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
