#' SYDNI Engine: sparse graphical synthesis
#'
#' @description
#' Fits subgroup-aware sparse Gaussian graphical models on feeder data and
#' generates synthetic draws with empirical marginal reconstruction. This
#' engine is intended when direct conditional-dependence structure is the
#' scientific target.
#'
#' @details
#' The graphical engine is SYDNI's sparse direct-dependence engine. It is
#' intended for settings where the analyst cares about conditional dependence,
#' network-like structure, or sparse relationships among synthesized variables.
#' Unlike the copula engine, which preserves dense joint association, the
#' graphical engine estimates a sparse latent precision structure.
#'
#' Let \eqn{D = \{x_i\}_{i=1}^{n_{\mathrm{obs}}}} be the observed data,
#' \eqn{V = \{v_1,\ldots,v_{p_V}\}} be the synthesized variables, and
#' \eqn{C = \{c_1,\ldots,c_q\}} be optional covariates. The requested synthetic
#' sample size is \eqn{n_{\mathrm{aug}}}. The engine supports numeric,
#' integer-like, factor, and ordered-factor synthesized columns.
#'
#' Covariates are handled through subgroup synthesis. If \code{covariates} are
#' supplied, each unique covariate combination defines a subgroup:
#' \deqn{
#'   g_i = \mathrm{interaction}(c_{i1},\ldots,c_{iq}).
#' }
#' If no covariates are supplied, all rows belong to a single group.
#'
#' The engine proceeds in eight steps:
#'
#' \enumerate{
#'   \item For each subgroup \eqn{g}, compute its observed size
#'     \eqn{n_g}, observed proportion \eqn{\pi_g = n_g / n_{\mathrm{obs}}},
#'     and feeder target size
#'     \deqn{
#'       n_g^{\mathrm{target}} =
#'       \max\left(
#'         n_g,\;
#'         \lceil \pi_g n_{\mathrm{aug}} \rceil,\;
#'         \max(10,\; 2p_V)
#'       \right).
#'     }
#'     The number of additional noisy feeder rows is
#'     \deqn{
#'       n_g^{\mathrm{syn}} =
#'       \max(n_g^{\mathrm{target}} - n_g,\; 0).
#'     }
#'
#'   \item Bootstrap \eqn{n_g^{\mathrm{syn}}} additional rows from each subgroup
#'     with replacement, inject noise into the synthesized columns only, and
#'     pool those noisy rows with the observed subgroup rows. Four noise methods
#'     are available through \code{noise_method}: \code{"discrete"},
#'     \code{"gaussian"}, \code{"truncnormal"}, and \code{"random"}. The
#'     method-specific distributions, automatic scaling rules for \code{c}, and
#'     type reconstruction details are implemented in
#'     \code{\link[=inject_noise]{inject_noise()}} and applied column-wise by
#'     \code{\link[=add_noise]{add_noise()}}. The pooled subgroup feeder data
#'     are denoted \eqn{D_g^{\mathrm{feed}}}.
#'
#'   \item Allocate final synthetic output rows across subgroups using the
#'     original subgroup proportions:
#'     \deqn{
#'       n_g^{\star} =
#'       \mathrm{round}(\pi_g n_{\mathrm{aug}}),
#'     }
#'     with any rounding discrepancy added to the largest allocated subgroup so
#'     that \eqn{\sum_g n_g^{\star} = n_{\mathrm{aug}}}.
#'
#'   \item Convert each subgroup feeder variable to a numeric structural score:
#'     factors and ordered factors are represented by integer level indices,
#'     and numeric or integer variables are represented numerically. These
#'     structural scores are transformed to a latent Gaussian scale. With
#'     \code{transform = "nonparanormal"}, the engine uses
#'     \code{huge::huge.npn()}; otherwise, it standardizes the structural matrix
#'     with \code{scale()}. Let \eqn{Z_g} denote the resulting latent matrix.
#'
#'   \item Identify active variables with finite non-zero latent standard
#'     deviation. Let \eqn{A_g \subseteq V} denote these active variables. The
#'     subgroup latent correlation matrix is initialized as an identity matrix
#'     over all synthesized variables, so inactive variables remain independent
#'     in the latent simulation unless active-variable estimates replace the
#'     corresponding block.
#'
#'   \item If \eqn{|A_g| \ge 2}, estimate a sparse latent precision matrix on
#'     \eqn{Z_g[A_g]} using graphical lasso through \pkg{huge}. The fitted
#'     precision matrix
#'     \deqn{
#'       \Omega_g
#'     }
#'     represents direct conditional dependence: zeros in \eqn{\Omega_g}
#'     correspond to absent edges in the active-variable latent Gaussian model.
#'
#'     If \code{selection = "fixed_lambda"}, the selected model is the fitted
#'     path value whose \eqn{\lambda} is closest to \code{lambda_fixed}. If
#'     \code{selection = "ebic"}, the model is selected with
#'     \code{huge::huge.select()} using an augmentation-adjusted EBIC penalty.
#'     Let
#'     \deqn{
#'       r_{\mathrm{aug}} =
#'       \frac{n_{\mathrm{aug}}}{n_{\mathrm{obs}}}.
#'     }
#'     The effective EBIC penalty is
#'     \deqn{
#'       \gamma_{\mathrm{eff}} =
#'       \min\left(
#'         1,\;
#'         \gamma_0 + 0.15 \log(1 + r_{\mathrm{aug}}) k
#'       \right),
#'     }
#'     where \eqn{\gamma_0} is \code{ebic_gamma}. Larger augmentation ratios
#'     therefore make graph selection more conservative.
#'
#'   \item Convert the selected precision matrix to a latent correlation matrix:
#'     \deqn{
#'       \Sigma_g^{A} =
#'       \mathrm{cov2cor}\left(\Omega_g^{-1}\right).
#'     }
#'     The active block is stabilized by symmetrizing it, replacing invalid
#'     entries, flooring eigenvalues below \code{pd_tol}, and rescaling to unit
#'     diagonal. If sparse fitting or selection fails, the engine falls back to
#'     the stabilized empirical latent correlation matrix for active variables.
#'     The active block is inserted into the full \eqn{p_V \times p_V} latent
#'     correlation matrix \eqn{\Sigma_g}, which is stabilized again before
#'     simulation.
#'
#'   \item Simulate latent Gaussian copula draws and reconstruct observed
#'     margins. For each subgroup,
#'     \deqn{
#'       Z_g^\star \sim \mathcal{N}(0,\; \Sigma_g),
#'       \qquad
#'       U_g = \Phi(Z_g^\star).
#'     }
#'     For variable \eqn{v}, synthetic values are reconstructed through the
#'     empirical quantile map of the subgroup feeder margin:
#'     \deqn{
#'       x_{iv}^{\star(g)} =
#'       \left(F_{vg}^{\mathrm{feed}}\right)^{-1}(u_{iv}),
#'     }
#'     where \eqn{F_{vg}^{\mathrm{feed}}} is the empirical distribution of
#'     variable \eqn{v} in \eqn{D_g^{\mathrm{feed}}}. Integer-like, factor, and
#'     ordered-factor variables are mapped back to valid observed support when
#'     \code{round_int = TRUE}.
#' }
#'
#' The final synthetic dataset is formed by stacking subgroup outputs. If
#' \code{covariates} are supplied, each synthetic subgroup block receives the
#' corresponding covariate values for that subgroup. The returned
#' \code{models} element stores, for each subgroup, the selected precision
#' matrix when available, the stabilized latent correlation matrix, the
#' \pkg{huge} selection object, the fitted \pkg{huge} path, fallback status, and
#' active variable names.
#'
#' This engine should be used when sparse direct dependence is substantively
#' meaningful; it is not a drop-in replacement for the dense
#' dependence-preserving copula engine.
#'
#' @param data A data frame or tibble.
#' @param variables Character vector of columns in `data` to synthesize.
#' @param covariates Character vector of columns in `data` to use as covariates.
#' @param noise_method Character. One of `"discrete"`, `"truncnormal"`,
#'   `"gaussian"`, or `"random"`.
#' @param c A numeric value in `(0, 1)` or `"auto"`.
#' @param round_int Logical. If `TRUE`, integer-like reconstructed values are
#'   rounded.
#' @param p Numeric. Probability control for discrete noise when `c = "auto"`.
#' @param f Numeric. Variance ratio control for continuous noise when `c =
#'   "auto"`.
#' @param r Numeric. Maximum noise magnitude for `noise_method = "random"`.
#' @param noise_bounds Optional numeric vector of length 2 used to clamp feeder
#'   rows before model fitting.
#' @param n_aug Numeric. Number of synthetic rows to generate.
#' @param k Numeric. SYDNI shrinkage exponent.
#' @param engine_control Optional named list. Supported elements include
#'   `transform`, `selection`, `ebic_gamma`, `lambda_fixed`, and `pd_tol`.
#' @param measurement_blocks Optional named list defining item blocks whose
#'   within-block latent correlations should be preserved while between-block
#'   latent correlations are shrunk. Also accepts a
#'   `cscore::composite_list()`-style list with `lower` and `higher` elements.
#' @param between_block_shrink Optional numeric value in `[0, 1]`. When
#'   `measurement_blocks` is supplied, between-block latent correlations are
#'   multiplied by this value. Defaults to `0.50`.
#' @param seed Numeric. Seed for reproducibility.
#'
#' @return A named list with elements `models` and `data`.
#'
#' @export
engine_graphical <- function(
    data,
    variables = NULL,
    covariates = NULL,
    noise_method = c("truncnormal", "gaussian", "discrete", "random"),
    c = "auto",
    round_int = TRUE,
    p = 0.5,
    f = 0.5,
    r = 1,
    noise_bounds = NULL,
    n_aug = nrow(data),
    k = 2,
    engine_control = NULL,
    measurement_blocks = NULL,
    between_block_shrink = NULL,
    seed = NA
) {
  noise_method <- match.arg(noise_method)
  engine_control <- validate_engine_control(engine_control)

  validate_engine_packages("graphical", "huge")

  if (is.null(variables)) {
    variables <- setdiff(colnames(data), covariates)
  }

  if (length(variables) < 1) {
    stop("`variables` must contain at least one synthesized column.")
  }

  validate_graphical_data(data, variables)

  control <- modifyList(
    list(
      transform = "nonparanormal",
      selection = "ebic",
      ebic_gamma = 0.5,
      lambda_fixed = NULL,
      pd_tol = 1e-8
    ),
    engine_control
  )

  normalized_blocks <- as_measurement_blocks(
    measurement_blocks = measurement_blocks,
    variables = variables
  )
  if (!is.null(normalized_blocks) && is.null(between_block_shrink)) {
    between_block_shrink <- 0.5
  }
  block_shrinkage_matrix <- if (!is.null(normalized_blocks)) {
    make_block_shrinkage_matrix(
      variables = variables,
      measurement_blocks = normalized_blocks,
      between_block_shrink = between_block_shrink
    )
  } else {
    NULL
  }

  group_plan <- build_group_plan(
    data = data,
    variables = variables,
    covariates = covariates,
    n_aug = n_aug
  )

  feed_data <- sample_group_feeder(
    group_plan = group_plan,
    variables = variables,
    noise_method = noise_method,
    c = c,
    p = p,
    f = f,
    r = r,
    round_int = round_int,
    noise_bounds = noise_bounds,
    seed = seed
  )

  group_n <- compute_group_output_sizes(group_plan$data, n_aug)
  aug_ratio <- n_aug / nrow(data)
  gamma_eff <- min(1, control$ebic_gamma + 0.15 * log1p(aug_ratio) * k)
  global_struct_mat <- make_structural_matrix(data, variables)
  global_latent_mat <- if (identical(control$transform, "nonparanormal")) {
    quiet_huge_npn(global_struct_mat)
  } else {
    scale(global_struct_mat)
  }
  global_latent_mat <- as.matrix(global_latent_mat)
  cor_global_reference <- stats::cor(global_latent_mat, use = "pairwise.complete.obs")
  cor_global_reference <- stabilize_correlation_matrix(cor_global_reference, tol = control$pd_tol)

  synth_predict <- stats::setNames(
    lapply(
      seq_along(group_plan$groups),
      function(i) {
        group_name <- group_plan$groups[[i]]
        data_g <- feed_data[feed_data$group_id == group_name, , drop = FALSE]
        n_sample <- as.integer(group_n[[group_name]])

        struct_mat <- make_structural_matrix(data_g, variables)
        latent_mat <- if (identical(control$transform, "nonparanormal")) {
          quiet_huge_npn(struct_mat)
        } else {
          scale(struct_mat)
        }
        latent_mat <- as.matrix(latent_mat)
        valid_sd <- apply(latent_mat, 2, stats::sd, na.rm = TRUE)
        active <- is.finite(valid_sd) & valid_sd > sqrt(.Machine$double.eps)

        sigma <- diag(length(variables))
        dimnames(sigma) <- list(variables, variables)
        graph_fit <- NULL
        select_fit <- NULL
        fallback <- FALSE

        if (sum(active) >= 2) {
          latent_active <- latent_mat[, active, drop = FALSE]

          graph_result <- tryCatch({
            fit <- huge::huge(
              x = latent_active,
              method = "glasso",
              verbose = FALSE
            )

            selected <- if (identical(control$selection, "fixed_lambda")) {
              if (is.null(control$lambda_fixed)) {
                stop("`engine_control$lambda_fixed` must be supplied when `selection = 'fixed_lambda'`.")
              }
              idx <- which.min(abs(fit$lambda - control$lambda_fixed))
              list(
                opt.icov = fit$icov[[idx]],
                opt.index = idx
              )
            } else {
              huge::huge.select(
                est = fit,
                criterion = "ebic",
                ebic.gamma = gamma_eff,
                verbose = FALSE
              )
            }

            omega_active <- selected$opt.icov
            sigma_active <- stats::cov2cor(solve(omega_active))
            sigma_active <- stabilize_correlation_matrix(sigma_active, tol = control$pd_tol)

            list(
              fit = fit,
              selected = selected,
              sigma = sigma_active,
              fallback = FALSE
            )
          }, error = function(e) {
            cor_active <- stats::cor(latent_active, use = "pairwise.complete.obs")
            cor_active <- stabilize_correlation_matrix(cor_active, tol = control$pd_tol)

            list(
              fit = NULL,
              selected = NULL,
              sigma = cor_active,
              fallback = TRUE,
              error = conditionMessage(e)
            )
          })

          sigma[active, active] <- graph_result$sigma
          graph_fit <- graph_result$fit
          select_fit <- graph_result$selected
          fallback <- isTRUE(graph_result$fallback)
        } else {
          fallback <- TRUE
        }
        sigma_raw <- stabilize_correlation_matrix(sigma, tol = control$pd_tol)
        sigma_adjusted <- sigma_raw
        if (!is.null(block_shrinkage_matrix)) {
          sigma_adjusted <- apply_block_shrinkage(sigma_adjusted, block_shrinkage_matrix)
          sigma_adjusted <- restore_within_block_correlations(
            cor_mat = sigma_adjusted,
            reference_cor = cor_global_reference,
            measurement_blocks = normalized_blocks
          )
          sigma_adjusted <- stabilize_correlation_matrix(sigma_adjusted, tol = control$pd_tol)
          sigma_adjusted <- restore_within_block_correlations(
            cor_mat = sigma_adjusted,
            reference_cor = cor_global_reference,
            measurement_blocks = normalized_blocks
          )
        }
        sigma <- sigma_adjusted
        sigma <- stabilize_correlation_matrix(sigma, tol = control$pd_tol)

        u_mat <- with_seed(seed + i - 1, simulate_gaussian_copula(n_sample, sigma))
        synth_df <- as.data.frame(u_mat)
        names(synth_df) <- variables

        for (j in seq_along(variables)) {
          synth_df[[variables[[j]]]] <- empirical_quantile_map(
            u = synth_df[[variables[[j]]]],
            x = data_g[[variables[[j]]]]
          )

          if (round_int) {
            synth_df[[variables[[j]]]] <- reconstruct_column(
              synth_df[[variables[[j]]]],
              data[[variables[[j]]]]
            )
          }
        }

        if (!is.null(covariates)) {
          cov_row <- group_plan$data[group_plan$data$group_id == group_name, covariates, drop = FALSE][1, , drop = FALSE]
          cov_block <- cov_row[rep(1, n_sample), , drop = FALSE]
          synth_df <- cbind(synth_df, cov_block)
        }

        list(
          model = list(
            precision = if (is.null(select_fit)) NULL else select_fit$opt.icov,
            correlation = sigma,
            correlation_raw = sigma_raw,
            correlation_measurement_adjusted = sigma_adjusted,
            correlation_global_reference = cor_global_reference,
            selection = select_fit,
            huge_fit = graph_fit,
            fallback = fallback,
            active_variables = variables[active],
            measurement_mode = if (!is.null(normalized_blocks)) "block_graphical" else NULL,
            measurement_blocks = normalized_blocks,
            between_block_shrink = between_block_shrink
          ),
          df = synth_df
        )
      }
    ),
    group_plan$groups
  )

  list(
    models = purrr::map(synth_predict, "model"),
    data = purrr::map_dfr(synth_predict, "df")
  )
}
