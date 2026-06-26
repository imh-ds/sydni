#' SYDNI Engine: Bayesian regression models
#'
#' @description
#' Fits a collection of Bayesian regression models with \pkg{brms} on the
#' observed data and draws posterior-predictive samples for a
#' bootstrap-resampled, noise-injected feeder dataset. Each target variable is
#' modeled conditionally on a user-chosen subset of other variables and optional
#' covariates.
#'
#' @details
#' The BRM engine is SYDNI's parametric conditional-model engine. It generates
#' synthetic data by fitting one Bayesian regression model per synthesized
#' variable on the observed data, then drawing posterior-predictive values for a
#' bootstrap-resampled and noise-injected feeder dataset. It is most appropriate
#' when the analyst wants explicit conditional relationships among variables and
#' wants the fitted models to remain inspectable through \pkg{brms}.
#'
#' Let \eqn{D = \{x_i\}_{i=1}^{n_{\mathrm{obs}}}} be the observed data,
#' \eqn{V = \{v_1,\ldots,v_p\}} be the synthesized variables, and
#' \eqn{C = \{c_1,\ldots,c_q\}} be optional covariates. The requested synthetic
#' sample size is \eqn{n_{\mathrm{aug}}}.
#'
#' The engine proceeds in six steps:
#'
#' \enumerate{
#'   \item Draw a feeder sample \eqn{\tilde{D}} of size
#'     \eqn{n_{\mathrm{aug}}} from \eqn{D} using \code{dplyr::slice_sample()}:
#'     \deqn{
#'       \tilde{D} = \{x_{s_i}: s_i \sim
#'       \mathrm{Uniform}(\{1,\ldots,n_{\mathrm{obs}}\})\}_{i=1}^{n_{\mathrm{aug}}}.
#'     }
#'     If \code{replace = FALSE}, the indices are sampled without replacement.
#'
#'   \item Inject noise into the synthesized columns of \eqn{\tilde{D}} to form
#'     the feeder data \eqn{\tilde{D}_{\varepsilon}}. Covariates are not
#'     perturbed. For a synthesized column \eqn{x}, the noisy value is
#'     \deqn{
#'       x_i^{\varepsilon} =
#'       \mathrm{clip}(x_i + \varepsilon_i,\; L_x,\; U_x),
#'     }
#'     where \eqn{[L_x,U_x]} is either \code{noise_bounds} or the observed range
#'     of \eqn{x}. Integer-like, factor, and ordered-factor columns are rounded
#'     and reconstructed after perturbation when \code{round_int = TRUE}.
#'     Four noise methods are available through \code{noise_method}:
#'     \code{"discrete"}, \code{"gaussian"}, \code{"truncnormal"}, and
#'     \code{"random"}. The method-specific distributions, automatic scaling
#'     rules for \code{c}, and type reconstruction details are implemented in
#'     \code{\link[=inject_noise]{inject_noise()}} and applied column-wise by
#'     \code{\link[=add_noise]{add_noise()}}.
#'
#'   \item For each target variable \eqn{v_m \in V}, choose a predictor set
#'     \eqn{Z_m \subseteq V \setminus \{v_m\}} using
#'     \code{predictor_method}. \code{"full"} uses all remaining synthesized
#'     variables. \code{"random"} samples \code{predictor_k} of them uniformly.
#'     \code{"information"} ranks candidates by normalized mutual information
#'     and takes the top \code{predictor_k}. \code{"random_information"} ranks
#'     candidates by normalized mutual information, keeps the top
#'     \eqn{5 k_{\mathrm{pred}}}, and samples \eqn{k_{\mathrm{pred}}} from
#'     that short list, where \eqn{k_{\mathrm{pred}}} is \code{predictor_k}.
#'
#'     The normalized mutual information score used for ranking is
#'     \deqn{
#'       \mathrm{NMI}(X,Y) =
#'       \frac{I(X;Y)}{\max\{H(X),H(Y)\}},
#'     }
#'     where \eqn{I(X;Y)} is mutual information and \eqn{H(\cdot)} is entropy.
#'     Numeric variables with many unique values are discretized into quantile
#'     bins before scoring. The final right-hand side for target \eqn{v_m} is
#'     \eqn{Z_m \cup C}; if this set is empty, the model is intercept-only:
#'     \deqn{
#'       v_m \sim 1.
#'     }
#'
#'   \item Assign a \pkg{brms} family \eqn{F_m} to each target. If
#'     \code{families} supplies a named family for \eqn{v_m}, that family is
#'     used. Otherwise, SYDNI chooses a family from the observed target type:
#'     binary factors use Bernoulli logit models, ordered factors and
#'     Likert-like integers use cumulative logit models, unordered factors use
#'     categorical logit models, non-negative integer counts use Poisson or
#'     negative-binomial models, proportions in \eqn{(0,1)} use beta models,
#'     and other numeric variables use Gaussian models.
#'
#'   \item Fit one conditional Bayesian regression model on the observed data
#'     for each target. In generic form, the linear predictor is
#'     \deqn{
#'       \eta_{im} = \alpha_m + z_{im}^{\top}\beta_m +
#'       c_i^{\top}\gamma_m,
#'     }
#'     with outcome model
#'     \deqn{
#'       v_{im} \mid z_{im}, c_i, \theta_m \sim F_m(\eta_{im}, \phi_m),
#'     }
#'     where \eqn{\theta_m} denotes all fitted parameters for target
#'     \eqn{v_m}; \eqn{\phi_m} denotes any family-specific auxiliary parameter
#'     such as residual scale, shape, or overdispersion.
#'
#'     Slope priors are shrunk as the requested augmentation grows:
#'     \deqn{
#'       \beta_{jm}, \gamma_{\ell m} \sim
#'       \mathcal{N}(0,\; \sigma_{\beta}^{2}),
#'       \qquad
#'       \sigma_{\beta} =
#'       \frac{\mathrm{base\_sd}}
#'       {\left(1 + n_{\mathrm{aug}} / n_{\mathrm{obs}}\right)^k}.
#'     }
#'     Larger values of \code{k} make the conditional models more conservative
#'     as the synthetic sample size grows. This shrinkage is intended to
#'     discourage synthetic augmentation from behaving like naive sample-size
#'     inflation.
#'
#'     The intercept prior is estimated from the observed target distribution:
#'     \deqn{
#'       \alpha_m \sim
#'       \mathcal{N}(\bar{v}_m,\; (s_m\,\mathrm{sd\_scale})^2),
#'     }
#'     where \eqn{\bar{v}_m} and \eqn{s_m} are the observed mean and standard
#'     deviation of target \eqn{v_m}. These priors are passed to \pkg{brms} as
#'     class \code{"b"} and class \code{"Intercept"} priors.
#'
#'   \item Generate synthetic target values from the posterior predictive
#'     distribution using the noisy feeder rows as \code{newdata}. For feeder
#'     row \eqn{i} and target \eqn{v_m},
#'     \deqn{
#'       v_{im}^{\star} \sim
#'       p(v_m \mid z_{im}^{\varepsilon}, c_i, D, F_m).
#'     }
#'     In implementation, \code{brms::posterior_predict()} returns a matrix
#'     \eqn{P_m} with posterior draws in rows and feeder observations in
#'     columns. SYDNI samples one posterior row independently for each feeder
#'     observation:
#'     \deqn{
#'       v_{im}^{\star} = P_m[r_i, i],
#'       \qquad
#'       r_i \sim \mathrm{Uniform}(\{1,\ldots,R_m\}),
#'     }
#'     where \eqn{R_m} is the number of retained posterior predictive draws.
#'     The resulting vector is reconstructed to match the observed target:
#'     factors are mapped back to valid levels, integer-like outcomes are
#'     rounded and clamped to the observed range, and continuous numeric
#'     outcomes are clamped to the observed range.
#'
#'     For Gaussian numeric targets, SYDNI instead draws from
#'     \code{brms::posterior_epred()} and adds centered empirical residuals.
#'     The residual donor row is shared across all targets when
#'     \code{measurement_blocks = NULL}; when measurement blocks are supplied,
#'     it is shared within each block and sampled independently across blocks.
#'     This preserves residual covariance within item batteries and avoids the
#'     reliability loss that can occur when each item receives independent
#'     posterior-predictive noise.
#' }
#'
#' The final synthetic dataset is assembled column-wise as
#' \deqn{
#'   D^{\star} =
#'   \left[v_1^{\star},\ldots,v_p^{\star}, C_{\mathrm{feeder}}\right],
#' }
#' where \eqn{C_{\mathrm{feeder}}} is included only when \code{covariates} are
#' supplied. The returned \code{models} element contains the fitted \pkg{brms}
#' model for each synthesized target, and the returned \code{data} element
#' contains \eqn{D^{\star}}.
#'
#' The BRM engine should be interpreted as a robustness and sensitivity-analysis
#' engine. It does not make synthetic rows into additional empirical evidence,
#' and it does not provide disclosure protection on its own.
#'
#' @param data A data frame or tibble.
#' @param variables Character vector of columns in `data` to synthesize.
#' @param covariates Character vector of columns in `data` to use as covariates.
#' @param noise_method Character. One of `"discrete"`, `"truncnormal"`,
#'   `"gaussian"`, or `"random"`.
#' @param c A numeric value in `(0, 1)` or `"auto"`.
#' @param round_int Logical. If `TRUE`, integer-like inputs are rounded after
#'   noise injection.
#' @param p Numeric. Probability control for discrete noise when `c = "auto"`.
#' @param f Numeric. Variance ratio control for continuous noise when `c =
#'   "auto"`.
#' @param r Numeric. Maximum noise magnitude for `noise_method = "random"`.
#' @param noise_bounds Optional numeric vector of length 2 used to clamp
#'   synthesized values.
#' @param predictor_method String. One of `"random_information"`, `"random"`,
#'   `"information"`, or `"full"`.
#' @param predictor_k Integer. Number of synthesized variables to use as
#'   predictors for each conditional model.
#' @param replace Logical. Determines whether bootstrap resampling uses
#'   replacement.
#' @param n_aug Numeric. Number of synthetic rows to generate.
#' @param families Optional named list of \pkg{brms} families to use manually.
#' @param base_sd Numeric. Base standard deviation for slope priors.
#' @param sd_scale Numeric. Scale multiplier used to estimate intercept priors.
#' @param k Numeric. Shrinkage exponent for slope priors.
#' @param chains Numeric. Number of Markov chains.
#' @param cores Integer or `"auto"`. Number of CPU cores for model fitting.
#' @param iter Numeric. Number of iterations per chain.
#' @param warmup Integer. Number of warmup iterations.
#' @param thin Integer. Thinning interval.
#' @param algorithm String. `brms` algorithm to use.
#' @param backend Character. One of `"rstan"` or `"cmdstanr"`.
#' @param measurement_blocks Optional named list defining item blocks whose
#'   within-block predictors should be preserved while between-block predictors
#'   are limited. Also accepts a `cscore::composite_list()`-style list with
#'   `lower` and `higher` elements.
#' @param between_block_shrink Optional numeric value in `[0, 1]`. When
#'   `measurement_blocks` is supplied, the between-block predictor count is
#'   `floor(predictor_k * between_block_shrink)` and between-block priors are
#'   tightened by the same multiplier. Defaults to `0.50`.
#' @param seed Numeric. Seed for reproducibility.
#'
#' @return A named list with elements `models` and `data`.
#'
#' @importFrom stats as.formula
#' @export
engine_brm <- function(
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
    families = NULL,
    base_sd = 10,
    sd_scale = 2,
    k = 2,
    chains = 2,
    cores = "auto",
    iter = 500,
    warmup = floor(iter / 2),
    thin = 2,
    algorithm = c("fullrank", "meanfield", "sampling", "pathfinder", "laplace", "fixed_param"),
    backend = c("rstan", "cmdstanr"),
    measurement_blocks = NULL,
    between_block_shrink = NULL,
    seed = NA
) {
  noise_method <- match.arg(noise_method)
  predictor_method <- match.arg(predictor_method)
  algorithm <- match.arg(algorithm)
  backend <- match.arg(backend)

  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("The 'brms' package must be installed to use `engine = 'brm'`.")
  }

  if (is.null(variables)) {
    variables <- setdiff(colnames(data), covariates)
  }

  if (length(variables) < 1) {
    stop("`variables` must contain at least one synthesized column.")
  }

  if (is.null(predictor_k)) {
    predictor_k <- max(1, round((length(variables) - 1) / 3))
    message(
      paste0("Predictor k not specified. Defaulting to ", predictor_k, ".")
    )
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

  variable_families <- assign_families(
    data = data,
    variables = variables,
    manual_families = families
  )

  formula_info <- with_seed(seed, stats::setNames(
    lapply(
      seq_along(variables),
      function(i) {
        target <- variables[[i]]
        predictor_roles <- if (is.null(normalized_blocks)) {
          predictors <- with_seed(
            seed + i - 1,
            select_predictors(
              target = target,
              data = data,
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
              data = data,
              variables = variables,
              predictor_method = predictor_method,
              predictor_k = predictor_k,
              measurement_blocks = normalized_blocks,
              between_block_shrink = between_block_shrink
            )
          )
        }
        predictors <- predictor_roles$predictors
        rhs_terms <- c(predictors, covariates)
        rhs <- if (length(rhs_terms) == 0) "1" else paste(rhs_terms, collapse = " + ")
        list(
          formula = stats::as.formula(paste0(target, " ~ ", rhs)),
          predictor_roles = predictor_roles
        )
      }
    ),
    variables
  ))
  formulas <- purrr::map(formula_info, "formula")

  pb <- progress::progress_bar$new(
    format = "  Generating Data [:bar] :percent ETA: :eta",
    total = length(variables),
    clear = FALSE,
    width = 100
  )

  b_shrinkage <- estimate_shrinkage_prior(
    n_obs = nrow(data),
    n_syn = n_aug,
    base_sd = base_sd,
    k = k
  )
  residual_donor_rows <- if (is.null(normalized_blocks)) {
    with_seed(
      seed + length(variables) + 1000L,
      sample.int(nrow(data), size = n_aug, replace = TRUE)
    )
  } else {
    stats::setNames(
      lapply(seq_along(normalized_blocks), function(i) {
        with_seed(
          seed + length(variables) + 1000L + i,
          sample.int(nrow(data), size = n_aug, replace = TRUE)
        )
      }),
      names(normalized_blocks)
    )
  }

  brm_models <- stats::setNames(
    lapply(
      seq_along(variables),
      function(i) {
        v <- variables[[i]]
        intercept_prior <- estimate_intercept_prior(
          y = data[[v]],
          sd_scale = sd_scale
        )

        if (isTRUE(getOption("sydni.brm_verbose", FALSE))) {
          message(
            paste0(
              "Fitting BRM target ", i, "/", length(variables), ": ", v,
              " | family = ", variable_families[[v]]$family,
              " | formula = ", paste(deparse(formulas[[v]]), collapse = "")
            )
          )
          utils::flush.console()
        }

        prior_info <- if (is.null(normalized_blocks)) {
          list(
            priors = b_shrinkage,
            meta = list(
              prior_mode = "global",
              same_block_sd = NA_real_,
              between_block_sd = NA_real_,
              covariate_sd = NA_real_,
              fallback_terms = character()
            )
          )
        } else {
          make_brm_role_priors(
            data = data,
            predictor_roles = formula_info[[v]]$predictor_roles,
            covariates = covariates,
            base_sd = base_sd,
            shrink_sd = base_sd / (1 + n_aug / nrow(data))^k,
            between_block_shrink = between_block_shrink
          )
        }

        model <- with_seed(seed + i - 1, suppressMessages(
          quiet_brm(
            formula = formulas[[v]],
            data = data,
            family = variable_families[[v]],
            prior = if (is.null(prior_info$priors)) intercept_prior else c(prior_info$priors, intercept_prior),
            chains = chains,
            cores = if (identical(cores, "auto")) parallel::detectCores() - 1 else cores,
            iter = iter,
            warmup = warmup,
            thin = thin,
            algorithm = algorithm,
            seed = seed,
            backend = backend,
            refresh = 0,
            silent = 2
          )
        ))

        preds <- feeder_data |>
          dplyr::select(
            dplyr::all_of(all.vars(formulas[[v]])[-1])
          )

        gaussian_numeric_target <- identical(variable_families[[v]]$family, "gaussian") &&
          (is.numeric(data[[v]]) || is.integer(data[[v]]))
        vec <- if (gaussian_numeric_target) {
          train_preds <- data |>
            dplyr::select(
              dplyr::all_of(all.vars(formulas[[v]])[-1])
            )
          epred_new <- brms::posterior_epred(model, newdata = preds)
          epred_train <- brms::posterior_epred(model, newdata = train_preds)
          donor_rows <- if (is.null(normalized_blocks)) {
            residual_donor_rows
          } else {
            residual_donor_rows[[formula_info[[v]]$predictor_roles$target_block]]
          }
          with_seed(seed + length(variables) + i - 1, sample_coherent_gaussian_epred_column(
            epred_new = epred_new,
            epred_train = epred_train,
            y_train = data[[v]],
            template = data[[v]],
            donor_rows = donor_rows
          ))
        } else {
          pp <- brms::posterior_predict(model, newdata = preds)
          with_seed(seed + length(variables) + i - 1, sample_posterior_column(
            pp = pp,
            template = data[[v]],
            family = variable_families[[v]]
          ))
        }

        pb$tick()

        list(
          model = model,
          value = vec,
          predictor_roles = list(
            within_block_predictors = formula_info[[v]]$predictor_roles$within_block_predictors,
            between_block_predictors = formula_info[[v]]$predictor_roles$between_block_predictors,
            covariates = covariates,
            target_block = formula_info[[v]]$predictor_roles$target_block,
            between_k = formula_info[[v]]$predictor_roles$between_k
          ),
          prior_roles = prior_info$meta
        )
      }
    ),
    variables
  )

  synthetic_data <- purrr::set_names(
    purrr::map_dfc(brm_models, "value"),
    variables
  )

  synth_models <- if (is.null(normalized_blocks)) {
    purrr::map(brm_models, "model")
  } else {
    purrr::map(brm_models, function(x) {
      list(
        model = x$model,
        predictor_roles = x$predictor_roles,
        prior_roles = x$prior_roles,
        measurement_mode = "block_predictor_brm",
        measurement_blocks = normalized_blocks,
        between_block_shrink = between_block_shrink
      )
    })
  }

  if (!is.null(covariates)) {
    synthetic_data <- cbind(
      synthetic_data,
      feeder_data[, covariates, drop = FALSE]
    )
  }

  list(
    models = synth_models,
    data = synthetic_data
  )
}
