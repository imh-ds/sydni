#' SYDNI Engine: latent factor synthesis
#'
#' @description
#' Fits subgroup-aware latent factor models on noisy feeder data and generates
#' synthetic psychometric indicators from the implied latent structure. This
#' engine is intended for ordinal or integer-like item batteries.
#'
#' @details
#' The factor engine is SYDNI's latent-structure engine for psychometric item
#' batteries and related ordinal or integer-like indicators. It is designed for
#' settings where the scientific target is shared latent measurement structure,
#' rather than only pairwise observed-variable association or target-wise
#' conditional prediction.
#'
#' Let \eqn{D = \{x_i\}_{i=1}^{n_{\mathrm{obs}}}} be the observed data,
#' \eqn{V = \{v_1,\ldots,v_{p_V}\}} be the synthesized item variables, and
#' \eqn{C = \{c_1,\ldots,c_q\}} be optional covariates. The requested synthetic
#' sample size is \eqn{n_{\mathrm{aug}}}. The engine first validates that every
#' synthesized column is factor, ordered factor, integer, or numeric
#' integer-like.
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
#'   \item Convert each subgroup feeder item to a numeric structural score:
#'     factors and ordered factors are represented by their integer level
#'     indices, and integer-like numeric items are represented numerically. An
#'     item is active in subgroup \eqn{g} if it has more than one distinct
#'     non-missing structural value. Let \eqn{A_g \subseteq V} denote the active
#'     items. Inactive items are later filled from their observed subgroup
#'     support.
#'
#'   \item If \eqn{|A_g| \ge 2}, estimate an item correlation matrix on active
#'     items. With \code{correlation_method = "polychoric"}, the engine first
#'     attempts
#'     \deqn{
#'       R_g = \mathrm{polychoric}(D_g^{\mathrm{feed}}[A_g]),
#'     }
#'     and falls back to pairwise Pearson correlation if polychoric estimation
#'     fails. Other correlation methods use pairwise Pearson correlation
#'     directly. The resulting matrix is stabilized before factor fitting.
#'
#'   \item Fit a factor model to \eqn{R_g} with \pkg{psych}. The effective
#'     number of factors is
#'     \deqn{
#'       q_g = \min(q_0,\; |A_g| - 1),
#'     }
#'     where \eqn{q_0} is \code{n_factors}, and extraction uses
#'     \code{extraction_method} with no rotation. For active subgroup items,
#'     the fitted factor model is represented as
#'     \deqn{
#'       Y_g = F_g \Lambda_g^\top + E_g,
#'     }
#'     where \eqn{F_g} is the latent factor score matrix, \eqn{\Lambda_g} is
#'     the loading matrix, and \eqn{E_g} contains item residuals with
#'     uniquenesses \eqn{\Psi_g}. If there are too few active items or factor
#'     fitting fails, the engine falls back to empirical subgroup
#'     reconstruction.
#'
#'   \item Regularize the fitted latent structure as augmentation grows. Let
#'     \deqn{
#'       r_{\mathrm{aug}} =
#'       \frac{n_{\mathrm{aug}}}{n_{\mathrm{obs}}}.
#'     }
#'     The loading multiplier is
#'     \deqn{
#'       m_\Lambda =
#'       \frac{1}
#'       {(1 + r_{\mathrm{aug}})^{k \cdot \mathrm{loading\_shrink}}},
#'     }
#'     so the shrunken loadings are
#'     \deqn{
#'       \Lambda_g^\star = m_\Lambda \Lambda_g.
#'     }
#'     If the fitted model supplies factor correlations \eqn{\Phi_g}, they are
#'     stabilized and regularized toward identity:
#'     \deqn{
#'       \Phi_g^\star =
#'       (1 - m_\Lambda) I + m_\Lambda \Phi_g.
#'     }
#'     If no factor-correlation matrix is supplied, \eqn{\Phi_g} is the identity
#'     matrix. Uniquenesses are lower-bounded by an augmentation-linked floor:
#'     \deqn{
#'       \psi_{jg}^\star =
#'       \min\left(
#'         0.99,\;
#'         \max\left[
#'           \psi_{jg},\;
#'           \min\{0.95,\; u_0(1 + \log(1 + r_{\mathrm{aug}}))\}
#'         \right]
#'       \right),
#'     }
#'     where \eqn{u_0} is \code{uniqueness_floor}.
#'
#'   \item Generate latent item values and reconstruct observed item scales. For
#'     each subgroup,
#'     \deqn{
#'       F_g^\star \sim \mathcal{N}(0,\; \Phi_g^\star),
#'       \qquad
#'       E_g^\star \sim \mathcal{N}(0,\; \Psi_g^\star),
#'     }
#'     then
#'     \deqn{
#'       Y_g^\star =
#'       F_g^\star(\Lambda_g^\star)^\top + E_g^\star.
#'     }
#'     Each latent item column is standardized, transformed to uniforms with
#'     the standard normal CDF, and mapped through the empirical subgroup feeder
#'     margin:
#'     \deqn{
#'       x_{ij}^{\star(g)} =
#'       \left(F_{jg}^{\mathrm{feed}}\right)^{-1}
#'       \left[
#'         \Phi\left(
#'           \frac{Y_{ij}^{\star(g)} - \bar{Y}_j^\star}
#'                {s_j^\star}
#'         \right)
#'       \right].
#'     }
#'     Integer-like values are reconstructed to match the observed storage type
#'     when \code{round_int = TRUE}. Inactive subgroup items are filled with the
#'     first observed non-missing subgroup value, then reconstructed to the
#'     observed type. If the fallback path is used, active items are generated
#'     by drawing \eqn{u_i \sim \mathrm{Uniform}(0,1)} and applying the same
#'     empirical quantile map.
#' }
#'
#' The final synthetic dataset is formed by stacking subgroup outputs. If
#' \code{covariates} are supplied, each synthetic subgroup block receives the
#' corresponding covariate values for that subgroup. The returned
#' \code{models} element stores, for each subgroup, loadings, factor
#' correlations, uniquenesses, the item correlation matrix, the \pkg{psych}
#' factor fit, fallback status, fallback reason, and active/inactive item sets.
#'
#' This engine is most appropriate for robustness checks involving Likert-style
#' item batteries, scales, and measurement models. It is not intended as a
#' general-purpose mixed-type generator.
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
#'   `n_factors`, `correlation_method`, `extraction_method`, `loading_shrink`,
#'   `uniqueness_floor`, and `threshold_method`.
#' @param seed Numeric. Seed for reproducibility.
#'
#' @return A named list with elements `models` and `data`.
#'
#' @export
engine_factor <- function(
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
    seed = NA
) {
  noise_method <- match.arg(noise_method)
  engine_control <- validate_engine_control(engine_control)

  validate_engine_packages("factor", "psych")

  if (is.null(variables)) {
    variables <- setdiff(colnames(data), covariates)
  }

  if (length(variables) < 1) {
    stop("`variables` must contain at least one synthesized column.")
  }

  validate_factor_data(data, variables)

  control <- modifyList(
    list(
      n_factors = 1L,
      correlation_method = "polychoric",
      extraction_method = "minres",
      loading_shrink = 1,
      uniqueness_floor = 0.1,
      threshold_method = "empirical"
    ),
    engine_control
  )

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
  loading_multiplier <- 1 / ((1 + aug_ratio) ^ (k * control$loading_shrink))
  uniqueness_floor <- min(0.95, control$uniqueness_floor * (1 + log1p(aug_ratio)))

  synth_predict <- stats::setNames(
    lapply(
      seq_along(group_plan$groups),
      function(i) {
        group_name <- group_plan$groups[[i]]
        data_g <- feed_data[feed_data$group_id == group_name, , drop = FALSE]
        n_sample <- as.integer(group_n[[group_name]])

        struct_df <- as.data.frame(lapply(data_g[, variables, drop = FALSE], coerce_structural_input))
        active <- vapply(
          struct_df,
          function(x) length(unique(stats::na.omit(x))) > 1,
          logical(1)
        )
        active_vars <- variables[active]
        inactive_vars <- variables[!active]

        synth_df <- as.data.frame(matrix(NA, nrow = n_sample, ncol = length(variables)))
        names(synth_df) <- variables

        fallback <- FALSE
        fallback_reason <- NULL
        loadings <- NULL
        phi <- NULL
        uniquenesses <- NULL
        cor_mat <- NULL
        fa_fit <- NULL

        if (length(inactive_vars) > 0) {
          warning(
            paste0(
              "The factor engine found zero-variance subgroup items in group '",
              group_name,
              "' and will reconstruct them from their observed subgroup support: ",
              paste(inactive_vars, collapse = ", "),
              "."
            )
          )
        }

        if (length(active_vars) >= 2) {
          struct_active <- struct_df[, active_vars, drop = FALSE]

          cor_mat <- if (identical(control$correlation_method, "polychoric")) {
            tryCatch(
              psych::polychoric(struct_active)$rho,
              error = function(e) stats::cor(struct_active, use = "pairwise.complete.obs")
            )
          } else {
            stats::cor(struct_active, use = "pairwise.complete.obs")
          }

          cor_mat <- stabilize_correlation_matrix(cor_mat)

          n_factors_eff <- min(as.integer(control$n_factors), length(active_vars) - 1L)

          fa_result <- tryCatch(
            psych::fa(
              r = cor_mat,
              nfactors = n_factors_eff,
              n.obs = nrow(struct_active),
              fm = control$extraction_method,
              rotate = "none",
              warnings = FALSE
            ),
            error = function(e) e
          )

          if (inherits(fa_result, "error")) {
            fallback <- TRUE
            fallback_reason <- conditionMessage(fa_result)
          } else {
            fa_fit <- fa_result

            loadings <- unclass(fa_fit$loadings)
            if (is.null(dim(loadings))) {
              loadings <- matrix(loadings, ncol = 1)
            }
            rownames(loadings) <- active_vars
            loadings <- loadings * loading_multiplier

            phi <- fa_fit$Phi
            if (is.null(phi)) {
              phi <- diag(ncol(loadings))
            }
            phi <- stabilize_correlation_matrix(phi)
            phi <- diag(1 - loading_multiplier, nrow = nrow(phi)) + phi * loading_multiplier
            phi <- stabilize_correlation_matrix(phi)

            uniquenesses <- fa_fit$uniquenesses
            uniquenesses <- pmin(0.99, pmax(uniquenesses, uniqueness_floor))
            names(uniquenesses) <- active_vars

            latent_f <- with_seed(seed + i - 1, simulate_latent_gaussian(n_sample, phi))
            residuals <- matrix(
              stats::rnorm(n_sample * nrow(loadings), sd = sqrt(uniquenesses)),
              nrow = n_sample,
              ncol = nrow(loadings),
              byrow = TRUE
            )
            latent_y <- latent_f %*% t(loadings) + residuals
            colnames(latent_y) <- active_vars

            for (v in active_vars) {
              latent_col <- latent_y[, v]
              latent_sd <- stats::sd(latent_col)
              z_col <- if (is.finite(latent_sd) && latent_sd > sqrt(.Machine$double.eps)) {
                (latent_col - mean(latent_col)) / latent_sd
              } else {
                rep(0, length(latent_col))
              }

              u <- stats::pnorm(z_col)
              synth_df[[v]] <- empirical_quantile_map(u = u, x = data_g[[v]])

              if (round_int) {
                synth_df[[v]] <- reconstruct_column(synth_df[[v]], data[[v]])
              }
            }
          }
        } else {
          fallback <- TRUE
          fallback_reason <- "Too few non-constant subgroup items remained for factor fitting."
        }

        if (fallback) {
          warning(
            paste0(
              "The factor engine fell back to empirical subgroup reconstruction for group '",
              group_name,
              "'. Reason: ",
              fallback_reason
            )
          )

          for (v in active_vars) {
            u <- with_seed(seed + 1000 + i, stats::runif(n_sample))
            synth_df[[v]] <- empirical_quantile_map(u = u, x = data_g[[v]])

            if (round_int) {
              synth_df[[v]] <- reconstruct_column(synth_df[[v]], data[[v]])
            }
          }
        }

        for (v in inactive_vars) {
          observed_nonmissing <- stats::na.omit(data_g[[v]])
          fill_value <- if (length(observed_nonmissing) == 0) NA else observed_nonmissing[[1]]
          synth_df[[v]] <- rep(fill_value, n_sample)

          if (round_int) {
            synth_df[[v]] <- reconstruct_column(synth_df[[v]], data[[v]])
          }
        }

        if (!is.null(covariates)) {
          cov_row <- group_plan$data[group_plan$data$group_id == group_name, covariates, drop = FALSE][1, , drop = FALSE]
          cov_block <- cov_row[rep(1, n_sample), , drop = FALSE]
          synth_df <- cbind(synth_df, cov_block)
        }

        list(
          model = list(
            loadings = loadings,
            phi = phi,
            uniquenesses = uniquenesses,
            correlation = cor_mat,
            factor_fit = fa_fit,
            fallback = fallback,
            fallback_reason = fallback_reason,
            active_variables = active_vars,
            inactive_variables = inactive_vars
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
