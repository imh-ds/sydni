#' SYDNI Engine: Bayesian additive regression trees
#'
#' @description
#' Fits a collection of regularized BART models on a bootstrap-resampled,
#' noise-injected feeder dataset and draws posterior predictions to create a
#' synthetic dataset. This engine is intended as the nonlinear conditional
#' synthesis path in SYDNI.
#'
#' @details
#' The BART engine is SYDNI's nonlinear conditional-model engine. It follows
#' the same broad feeder-data logic as \code{engine_brm()}, but replaces
#' parametric Bayesian regressions with regularized Bayesian additive regression
#' trees from \pkg{dbarts}. It is intended for settings where nonlinearities or
#' interactions among predictors may be important.
#'
#' Let \eqn{D = \{x_i\}_{i=1}^{n_{\mathrm{obs}}}} be the observed data,
#' \eqn{V = \{v_1,\ldots,v_p\}} be the synthesized variables, and
#' \eqn{C = \{c_1,\ldots,c_q\}} be optional covariates. The requested synthetic
#' sample size is \eqn{n_{\mathrm{aug}}}.
#'
#' The engine proceeds in seven steps:
#'
#' \enumerate{
#'   \item Draw a feeder sample \eqn{\tilde{D}} of size
#'     \eqn{n_{\mathrm{aug}}} from \eqn{D} using \code{dplyr::slice_sample()}.
#'     Sampling uses \code{replace} exactly as supplied.
#'
#'   \item Inject noise into the synthesized columns of \eqn{\tilde{D}} to form
#'     feeder data \eqn{\tilde{D}_{\varepsilon}}. Covariates are not perturbed.
#'     Four noise methods are available through \code{noise_method}:
#'     \code{"discrete"}, \code{"gaussian"}, \code{"truncnormal"}, and
#'     \code{"random"}. The method-specific distributions, automatic scaling
#'     rules for \code{c}, and type reconstruction details are implemented in
#'     \code{\link[=inject_noise]{inject_noise()}} and applied column-wise by
#'     \code{\link[=add_noise]{add_noise()}}.
#'
#'   \item Classify each target from its observed storage type and support.
#'     Binary targets are detected as two-level factors or two-valued vectors.
#'     Factors, ordered factors, integers, and numeric columns whose non-missing
#'     values are all effectively whole numbers are treated as ordinal-like
#'     numeric approximations during fitting. All other numeric columns are
#'     treated as continuous.
#'
#'   \item For each target variable \eqn{v_m \in V}, choose a predictor set
#'     \eqn{Z_m \subseteq V \setminus \{v_m\}} from the feeder data using
#'     \code{predictor_method}. \code{"full"} uses all remaining synthesized
#'     variables. \code{"random"} samples \code{predictor_k} of them uniformly.
#'     \code{"information"} ranks candidates by normalized mutual information
#'     and takes the top \code{predictor_k}. \code{"random_information"} ranks
#'     candidates by normalized mutual information, keeps the top
#'     \eqn{5 k_{\mathrm{pred}}}, and samples \eqn{k_{\mathrm{pred}}} from
#'     that short list, where \eqn{k_{\mathrm{pred}}} is \code{predictor_k}.
#'     Covariates are then appended to the selected predictor set.
#'
#'     The normalized mutual information score used for ranking is
#'     \deqn{
#'       \mathrm{NMI}(X,Y) =
#'       \frac{I(X;Y)}{\max\{H(X),H(Y)\}},
#'     }
#'     where \eqn{I(X;Y)} is mutual information and \eqn{H(\cdot)} is entropy.
#'     Numeric variables with many unique values are discretized into quantile
#'     bins before scoring.
#'
#'   \item Convert the selected predictors and covariates into a numeric design
#'     matrix using \code{model.matrix(~ . - 1)}. If no predictors or covariates
#'     are available, an intercept-only frame is used. For feeder row \eqn{i},
#'     let \eqn{x_{im}^{\varepsilon}} denote this numeric predictor vector for
#'     target \eqn{v_m}.
#'
#'   \item Fit one BART model per target on the noisy feeder data. In schematic
#'     form, BART represents the conditional mean as a sum of regression trees:
#'     \deqn{
#'       y_{im}^{\varepsilon} =
#'       \sum_{t = 1}^{T} g_{tm}(x_{im}^{\varepsilon}) +
#'       e_{im},
#'     }
#'     where \eqn{T} is \code{num_trees}, each \eqn{g_{tm}} is a regression
#'     tree, and \eqn{e_{im}} is the residual term implied by the BART fit.
#'     Factor and ordered targets are converted to integer scores before
#'     fitting; binary targets are also fit numerically.
#'
#'     SYDNI makes the BART prior more conservative as the augmentation ratio
#'     grows. Let
#'     \deqn{
#'       r_{\mathrm{aug}} =
#'       \frac{n_{\mathrm{aug}}}{n_{\mathrm{obs}}}.
#'     }
#'     The effective BART regularization parameters passed to
#'     \code{dbarts::bart()} are
#'     \deqn{
#'       k_{\mathrm{eff}} =
#'       k_{\mathrm{bart}} (1 + r_{\mathrm{aug}})^k,
#'     }
#'     \deqn{
#'       \mathrm{power}_{\mathrm{eff}} =
#'       \mathrm{power} + \log(1 + r_{\mathrm{aug}}),
#'     }
#'     and
#'     \deqn{
#'       \mathrm{base}_{\mathrm{eff}} =
#'       \min(0.99,\; \mathrm{base} + 0.02 \log(1 + r_{\mathrm{aug}})).
#'     }
#'     Here \code{k_bart}, \code{power}, and \code{base} come from
#'     \code{engine_control}. Increasing \eqn{k_{\mathrm{eff}}} shrinks
#'     terminal-node effects, while larger \eqn{\mathrm{power}_{\mathrm{eff}}}
#'     and \eqn{\mathrm{base}_{\mathrm{eff}}} make deeper trees less likely.
#'
#'   \item Draw posterior predictions for the feeder design matrix. The
#'     \pkg{dbarts} fit returns a matrix \eqn{P_m} of test-set posterior
#'     predictions with posterior draws in rows and feeder observations in
#'     columns. SYDNI samples one posterior draw independently for each feeder
#'     observation:
#'     \deqn{
#'       \tilde{v}_{im} = P_m[r_i, i],
#'       \qquad
#'       r_i \sim \mathrm{Uniform}(\{1,\ldots,R_m\}).
#'     }
#'     By default, these conditional-mean draws are converted to posterior
#'     predictive draws by adding centered empirical residuals sampled from
#'     \eqn{y_m - \bar{P}_m}. When no measurement blocks are supplied, the same
#'     residual-donor row is used across targets for a synthetic observation, so
#'     residual variation preserves cross-item covariance instead of adding
#'     independent item-wise noise. When measurement blocks are supplied, the
#'     donor row is shared within each block but sampled independently across
#'     blocks, preserving measurement families while allowing between-block
#'     relations to remain conservative. This step is important for ordinal item
#'     batteries: without residual variation, BART can return conditional means
#'     clustered near the middle of the response scale, and rounding can collapse
#'     a synthesized item to a single category. With independent residual
#'     variation, item reliabilities can be artificially damaged. The amount of
#'     residual variation is controlled by
#'     \code{engine_control$residual_scale}. Set
#'     \code{engine_control$posterior_predictive = FALSE} to recover the
#'     conditional-mean-only behavior.
#'
#'     For binary targets, posterior-predictive mode treats the sampled BART
#'     value as a probability after clamping to \eqn{[0,1]} and draws a
#'     Bernoulli outcome. If posterior-predictive mode is disabled, the sampled
#'     value is thresholded as
#'     \deqn{
#'       v_{im}^{\star} =
#'       \mathbb{I}\{\tilde{v}_{im} \ge \tau\},
#'     }
#'     where \eqn{\tau} is \code{binary_threshold}. All targets are then
#'     reconstructed to match the observed support: factors are mapped back to
#'     valid levels, integer-like outcomes are rounded and clamped to the
#'     observed range, and continuous numeric outcomes are clamped to the
#'     observed range.
#' }
#'
#' The final synthetic dataset is assembled column-wise across target
#' variables. If \code{covariates} are supplied, the covariate columns from the
#' feeder data are appended unchanged. The returned \code{models} element
#' contains the fitted \pkg{dbarts} model for each synthesized target.
#'
#' This engine is useful for nonlinear conditional robustness checks, but it is
#' not a general privacy mechanism and should not be used to treat synthetic
#' rows as additional empirical evidence.
#'
#' @param data A data frame or tibble.
#' @param variables Character vector of columns in `data` to synthesize.
#' @param covariates Character vector of columns in `data` to use as covariates.
#' @param noise_method Character. One of `"discrete"`, `"truncnormal"`,
#'   `"gaussian"`, or `"random"`.
#' @param c A numeric value in `(0, 1)` or `"auto"`.
#' @param round_int Logical. If `TRUE`, integer-like inputs are rounded after
#'   reconstruction.
#' @param p Numeric. Probability control for discrete noise when `c = "auto"`.
#' @param f Numeric. Variance ratio control for continuous noise when `c =
#'   "auto"`.
#' @param r Numeric. Maximum noise magnitude for `noise_method = "random"`.
#' @param noise_bounds Optional numeric vector of length 2 used to clamp feeder
#'   rows before model fitting.
#' @param predictor_method String. One of `"random_information"`, `"random"`,
#'   `"information"`, or `"full"`.
#' @param predictor_k Integer. Number of synthesized variables to use as
#'   predictors for each conditional model.
#' @param replace Logical. Determines whether bootstrap resampling uses
#'   replacement.
#' @param n_aug Numeric. Number of synthetic rows to generate.
#' @param k Numeric. SYDNI shrinkage exponent. Used to scale BART regularization
#'   as augmentation grows.
#' @param engine_control Optional named list. Supported elements include
#'   `num_trees`, `num_burn`, `num_draws`, `k_bart`, `power`, `base`,
#'   `binary_threshold`, `posterior_predictive`, and `residual_scale`.
#' @param measurement_blocks Optional named list defining item blocks whose
#'   within-block predictors should be preserved while between-block predictors
#'   are limited. Also accepts a `cscore::composite_list()`-style list with
#'   `lower` and `higher` elements.
#' @param between_block_shrink Optional numeric value in `[0, 1]`. When
#'   `measurement_blocks` is supplied, the between-block predictor count is
#'   `floor(predictor_k * between_block_shrink)`. Defaults to `0.50`.
#' @param seed Numeric. Seed for reproducibility.
#'
#' @return A named list with elements `models` and `data`.
#'
#' @importFrom utils modifyList
#' @export
engine_bart <- function(
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
    predictor_method = c("random_information", "random", "information", "full"),
    predictor_k = NULL,
    replace = TRUE,
    n_aug = nrow(data),
    k = 2,
    engine_control = NULL,
    measurement_blocks = NULL,
    between_block_shrink = NULL,
    seed = NA
) {
  noise_method <- match.arg(noise_method)
  predictor_method <- match.arg(predictor_method)
  engine_control <- validate_engine_control(engine_control)

  validate_engine_packages("bart", "dbarts")

  if (is.null(variables)) {
    variables <- setdiff(colnames(data), covariates)
  }

  if (length(variables) < 1) {
    stop("`variables` must contain at least one synthesized column.")
  }

  bart_targets <- classify_bart_targets(data, variables)

  if (length(bart_targets$ordinal_approx) > 0) {
    warning(
      paste0(
        "`bart` identified ordinal-like synthesized columns and will treat them as numeric during fitting, then round/clamp them back to the observed support: ",
        paste(bart_targets$ordinal_approx, collapse = ", "),
        "."
      )
    )
  }

  if (length(bart_targets$numeric_ordinal) > 0) {
    warning(
      paste0(
        "The following synthesized columns are not stored as factor/ordered class but look ordinal-like, so `bart` is treating them as numeric approximations: ",
        paste(bart_targets$numeric_ordinal, collapse = ", "),
        "."
      )
    )
  }

  if (is.null(predictor_k)) {
    predictor_k <- max(1, round((length(variables) - 1) / 3))
  }

  predictor_k <- min(predictor_k, max(length(variables) - 1, 0))
  normalized_blocks <- as_measurement_blocks(
    measurement_blocks = measurement_blocks,
    variables = variables
  )
  if (!is.null(normalized_blocks) && is.null(between_block_shrink)) {
    between_block_shrink <- 0.5
  }
  if (!is.null(normalized_blocks)) {
    make_block_shrinkage_matrix(
      variables = variables,
      measurement_blocks = normalized_blocks,
      between_block_shrink = between_block_shrink
    )
  }

  control <- modifyList(
    list(
      num_trees = 50L,
      num_burn = 250L,
      num_draws = 250L,
      k_bart = 2,
      power = 2,
      base = 0.95,
      binary_threshold = 0.5,
      posterior_predictive = TRUE,
      residual_scale = 1
    ),
    engine_control
  )
  if (
    !is.numeric(control$residual_scale) ||
      length(control$residual_scale) != 1 ||
      is.na(control$residual_scale) ||
      control$residual_scale < 0
  ) {
    stop("`engine_control$residual_scale` must be a single non-negative numeric value.")
  }

  sampled_data <- with_seed(seed, dplyr::slice_sample(
    .data = data,
    replace = replace,
    n = n_aug
  ))

  feeder_data <- add_noise(
    data = sampled_data,
    variables = dplyr::all_of(variables),
    method = noise_method,
    c = c,
    p = p,
    f = f,
    round_int = round_int,
    r = r,
    bounds = noise_bounds
  )

  aug_ratio <- n_aug / nrow(data)
  effective_k <- control$k_bart * ((1 + aug_ratio) ^ k)
  effective_power <- control$power + log1p(aug_ratio)
  effective_base <- min(0.99, control$base + 0.02 * log1p(aug_ratio))
  residual_donor_rows <- if (is.null(normalized_blocks)) {
    with_seed(
      seed + length(variables) + 1000L,
      sample.int(nrow(feeder_data), size = n_aug, replace = TRUE)
    )
  } else {
    stats::setNames(
      lapply(seq_along(normalized_blocks), function(i) {
        with_seed(
          seed + length(variables) + 1000L + i,
          sample.int(nrow(feeder_data), size = n_aug, replace = TRUE)
        )
      }),
      names(normalized_blocks)
    )
  }

  bart_models <- stats::setNames(
    lapply(
      seq_along(variables),
      function(i) {
        target <- variables[[i]]
        predictor_roles <- if (is.null(normalized_blocks)) {
          predictors <- with_seed(
            seed + i - 1,
            select_predictors(
              target = target,
              data = feeder_data,
              variables = variables,
              predictor_method = predictor_method,
              predictor_k = predictor_k
            )
          )
          list(
            predictors = predictors,
            within_block_predictors = character(),
            between_block_predictors = predictors,
            between_k = predictor_k,
            target_block = NA_character_
          )
        } else {
          with_seed(
            seed + i - 1,
            select_measurement_predictors(
              target = target,
              data = feeder_data,
              variables = variables,
              predictor_method = predictor_method,
              predictor_k = predictor_k,
              measurement_blocks = normalized_blocks,
              between_block_shrink = between_block_shrink
            )
          )
        }
        predictors <- predictor_roles$predictors

        predictor_cols <- c(predictors, covariates)
        x_frame <- make_predictor_frame(feeder_data, predictor_cols)
        x_train <- make_numeric_design_matrix(x_frame)
        x_test <- x_train

        y_train <- feeder_data[[target]]
        target_kind <- bart_targets$kinds[[target]]
        binary_target <- identical(target_kind, "binary")
        y_numeric <- if (is.factor(y_train) || is.ordered(y_train)) {
          as.integer(y_train)
        } else {
          as.numeric(y_train)
        }

        fit <- with_seed(seed + i - 1, dbarts::bart(
          x.train = x_train,
          y.train = y_numeric,
          x.test = x_test,
          ntree = as.integer(control$num_trees),
          ndpost = as.integer(control$num_draws),
          nskip = as.integer(control$num_burn),
          keepevery = 1L,
          keeptrainfits = FALSE,
          verbose = FALSE,
          k = effective_k,
          power = effective_power,
          base = effective_base
        ))

        pp <- fit$yhat.test
        if (is.null(pp)) {
          stop("`dbarts::bart()` did not return test-set posterior draws.")
        }

        vec <- with_seed(seed + length(variables) + i - 1, {
          rows <- sample.int(nrow(pp), ncol(pp), replace = TRUE)
          raw <- pp[cbind(rows, seq_len(ncol(pp)))]

          if (binary_target) {
            if (isTRUE(control$posterior_predictive)) {
              prob <- pmin(pmax(raw, 0), 1)
              raw <- stats::rbinom(length(prob), size = 1, prob = prob)
            } else {
              raw <- as.numeric(raw >= control$binary_threshold)
            }
          } else if (isTRUE(control$posterior_predictive)) {
            fitted_mean <- colMeans(pp, na.rm = TRUE)
            residuals <- y_numeric - fitted_mean
            residuals[!is.finite(residuals)] <- 0
            if (length(residuals) > 0 && stats::sd(residuals) > sqrt(.Machine$double.eps)) {
              residuals <- residuals - mean(residuals)
              donor_rows <- if (is.null(normalized_blocks)) {
                residual_donor_rows
              } else {
                residual_donor_rows[[predictor_roles$target_block]]
              }
              raw <- raw + residuals[donor_rows] * control$residual_scale
            }
          }

          reconstruct_column(raw, data[[target]])
        })

        list(
          model = fit,
          value = vec,
          predictors = predictor_cols,
          predictor_roles = list(
            within_block_predictors = predictor_roles$within_block_predictors,
            between_block_predictors = predictor_roles$between_block_predictors,
            covariates = covariates,
            target_block = predictor_roles$target_block,
            between_k = predictor_roles$between_k
          )
        )
      }
    ),
    variables
  )

  synthetic_data <- purrr::set_names(
    purrr::map_dfc(bart_models, "value"),
    variables
  )

  if (!is.null(covariates)) {
    synthetic_data <- cbind(
      synthetic_data,
      feeder_data[, covariates, drop = FALSE]
    )
  }

  list(
    models = if (is.null(normalized_blocks)) {
      purrr::map(bart_models, "model")
    } else {
      purrr::map(bart_models, function(x) {
        list(
          model = x$model,
          predictors = x$predictors,
          predictor_roles = x$predictor_roles,
          measurement_mode = "block_predictor_bart",
          measurement_blocks = normalized_blocks,
          between_block_shrink = between_block_shrink
        )
      })
    },
    data = synthetic_data
  )
}
