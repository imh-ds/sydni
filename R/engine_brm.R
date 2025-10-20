#' SYDNI Engine: Bayesian Regression Model
#'
#' @description Fits a collection of Bayesian regression models with the
#'   \code{brms} package and draws posterior-predictive samples to create a
#'   fully synthetic data set.  Each target variable is modeled conditionally on
#'   a user-chosen subset of other variables and (optionally) external
#'   covariates, allowing realistic multivariate dependence. The function
#'   supports several noise-injection schemes, adaptive shrinkage priors, and
#'   multiple Stan backends as a general engine for the SYDNI framework.
#'
#'
#' @details
#'
#' \code{engine_brm()} is a two-stage synthesizer:
#'
#' \strong{Perturb the original data.} A bootstrap sample of size
#' \eqn{n_{\mathrm{aug}}} is drawn from the input data (with or without
#' replacement). Before any modelling, the selected columns are jittered by an
#' i.i.d. noise term \eqn{\varepsilon}. This step both widens the support of
#' each variable (to allow generative flexibility) and introduces uncertainty.
#'
#' The noise distribution depends on the variable type (continuous vs.
#' integer-like) and the \code{noise_method}:
#'
#'     \itemize{
#'       \item \code{"discrete"}:
#'         \emph{If} the variable is integer-like, adds noise from
#'         \eqn{\varepsilon \in \{-1, 0, +1\}} with
#'         \eqn{P(\varepsilon = \pm 1) = c / 2}, \eqn{P(\varepsilon = 0) = 1 - c}.
#'         \emph{If} the variable is continuous, instead adds
#'         \eqn{\varepsilon \sim \mathrm{Uniform}(-c \cdot \mathrm{sd}(x),\ c \cdot \mathrm{sd}(x))}.
#'
#'       \item \code{"gaussian"}:
#'         Adds Gaussian noise
#'         \eqn{\varepsilon \sim \mathcal{N}(0,\ (c \cdot \mathrm{sd}(x))^2)}.
#'         Rounded if the input is integer-like.
#'
#'       \item \code{"truncnormal"}:
#'         Adds truncated normal noise
#'         \eqn{\varepsilon \sim \mathrm{TruncNormal}(0,\ (c \cdot \mathrm{sd}(x))^2,\ [-c \cdot \mathrm{sd}(x),\ c \cdot \mathrm{sd}(x)])}.
#'         Also rounded if the input is integer-like.
#'
#'       \item \code{"random"}:
#'         Adds unscaled uniform noise
#'         \eqn{\varepsilon \sim \mathrm{Uniform}(-r,\ r)},
#'         independent of the data's standard deviation.
#'     }
#'
#' For all methods, the post-noise value is clamped to either the original or
#' user-specified bounds, and optionally rounded (via \code{round_int}) to
#' preserve the data's structure.
#'
#' \strong{Model each target column conditionally on others.} For each variable
#' \eqn{y_m} listed in \code{variables}, a subset of predictors \eqn{z_m} is
#' selected:
#'
#'     \itemize{
#'       \item \code{"full"}: all other variables
#'       \item \code{"random"}: a random subset of size \code{predictor_k}
#'       \item \code{"information"}: top \eqn{k} variables ranked by mutual information with \eqn{y_m}
#'       \item \code{"random_information"}: random draw from the top \eqn{5k} informative variables
#'     }
#'
#' Any covariates passed via \code{covariates} are also included in the model.
#'
#' Each \eqn{y_m} is modelled with a Bayesian regression of the form:
#'
#'     \deqn{
#'       \eta_m = \beta_{0m} + z_m^T \beta_m
#'     }
#'
#' The outcome distribution is chosen by \code{brms} based on variable type (or
#' manually via \code{families}).
#'
#' \strong{Adaptive shrinkage prior:}  
#' Each regression coefficient \eqn{\beta_{jm}} is given a normal prior whose
#' standard deviation shrinks based on the relative size of the synthetic sample.
#' Specifically, the prior is defined as:
#'
#' \deqn{
#'   \beta_{jm} \sim \mathcal{N}(0,\ \sigma_\beta^2), \quad
#'   \sigma_\beta = \frac{\mathrm{base\_sd}}{(1 + \frac{n_{\mathrm{syn}}}{n_{\mathrm{obs}}})^k}
#' }
#'
#' where \eqn{n_{\mathrm{obs}}} is the number of original (observed) data points,
#' \eqn{n_{\mathrm{syn}}} is the number of synthetic records, and \eqn{k} is a
#' scaling exponent controlling how aggressively the prior tightens.
#'
#' As \eqn{n_{\mathrm{syn}}} increases relative to \eqn{n_{\mathrm{obs}}},
#' the prior variance shrinks, discouraging the model from overfitting
#' synthetic noise. This helps maintain stability during iterative generation
#' or large-scale synthesis.
#'
#' \strong{Posterior prediction:} From each fitted model \eqn{\mathcal{M}_m}, a
#' posterior-predictive draw is made for each row in the noisy data:
#'
#'     \deqn{
#'       \hat{y}_m^{(i)} \sim p(y_m \mid z_m^{(i)},\ \mathcal{M}_m), \quad i = 1, \dots, n_{\mathrm{aug}}
#'     }
#'
#' These predicted values constitute the synthetic version of \eqn{y_m}.
#'
#' @param data A data frame or tibble.
#' @param variables A tidyselect specification indicating which columns in
#'   \code{data} to synthesize.
#' @param covariates A tidyselect specification indicating which columns in
#'   \code{data} to use as covariates in the synthetic data generation process.
#' @param noise_method Character. One of "discrete", "truncnormal", "gaussian",
#'   or "random".
#'   \describe{
#'     \item{\code{discrete}}{
#'       For integers or ordered factors, adds noise with possible values
#'       -1, 0, or +1, where
#'       \eqn{P(\pm1) = \frac{c}{2}, \quad P(0) = 1 - c}.
#'       For continuous inputs, adds uniform noise:
#'       \eqn{\varepsilon \sim \mathrm{Uniform}(-c \cdot \mathrm{sd}(x),\; c \cdot \mathrm{sd}(x))}.
#'     }
#'     \item{\code{truncnormal}}{
#'       Adds noise from a truncated normal distribution centered at 0:
#'       \eqn{\varepsilon \sim \mathrm{TruncNormal}(0, (c \cdot \mathrm{sd}(x))^2,\;
#'       [-c \cdot \mathrm{sd}(x),\; c \cdot \mathrm{sd}(x)])}.
#'       Rounded to ±1 if input is integer-like.
#'     }
#'     \item{\code{gaussian}}{
#'       Adds noise from a normal distribution:
#'       \eqn{\varepsilon \sim \mathcal{N}(0,\; (c \cdot \mathrm{sd}(x))^2)}.
#'       Rounded if input is integer-like.
#'     }
#'     \item{\code{random}}{
#'       Adds uniform noise independent of data scale:
#'       \eqn{\varepsilon \sim \mathrm{Uniform}(-r,\; r)}.
#'     }
#'   }
#' @param c A numeric value in (0, 1). Used to scale the standard deviation of
#'   the noise. Defaults to "auto", in which case it is computed based on either
#'   \code{p} (for discrete) or \code{f} (for continuous) inputs.
#' @param round_int Logical. If TRUE (default), integer or ordered factor inputs
#'   are rounded to the nearest whole number after noise is added.
#' @param p Numeric. Probability mass to assign to jumps of ±1 in the discrete
#'   case (used only if \code{c = "auto"}).
#' @param f Numeric. Fraction of standard deviation to approximate when using
#'   continuous noise (used only if \code{c = "auto"}).
#' @param r Numeric. Maximum magnitude of noise when using \code{method =
#'   "random"}.
#' @param noise_bounds Optional numeric vector of length 2. If specified, clamps
#'   output to these lower and upper bounds. Otherwise, clamps to the range of
#'   \code{x}.
#' @param predictor_method String. Specifies the method used to determine which
#'   variables to use in the synthesis engine as predictors. Default is
#'   \code{random_information}.
#' @param predictor_k Integer. Specifies the number of variables to use in the
#'   synthesis engine as predictors. Default is \code{NULL}, which automatically
#'   determines the number.
#' @param replace Logic. Determines whether to resample with replacement.
#'   Default is \code{TRUE}.
#' @param n_aug Numeric. The sample size of the synthetic data to generate.
#' @param families A named list of families to use manually, where names
#'   correspond to column names in \code{data} and values are \code{brmsfamily}
#'   objects.
#' @param base_sd Numeric. The base standard deviation for the prior when no
#'   synthetic data is used. Default is 10.
#' @param sd_scale Numeric. Multiplier for the standard deviation of `y` to set
#'   the prior's spread. A value of 2 approximates a 95% coverage around the
#'   sample mean. Default is 2.
#' @param k Numeric. A scaling exponent controlling the rate of shrinkage.
#'   Higher values shrink the prior more aggressively as synthetic data
#'   increases. Default is 2.
#' @param chains Numeric. Number of Markov chains (Default is 2).
#' @param cores Integer. Number of CPU cores to use for parallel computation. By
#'   default, the number of available cores minus one is used:
#'   \code{parallel::detectCores() - 1}. For best performance, you may set this
#'   manually to match your system's capabilities.
#' @param iter Numeric. Number of iterations per chain, inclusive of warmup.
#'   Default is 500.
#' @param warmup Integer. Specifies the number of warmup (or burn-in)
#'   iterations. These iterations are used for step size adaptation and are not
#'   included in the posterior samples. By default, this is set to
#'   \code{floor(iter / 2)}.
#' @param thin Integer. Specifies the thinning interval for saving iterations.
#'   Every \code{thin}-th iteration is retained. Thinning helps reduce memory
#'   usage and autocorrelation in long chains. Defaults to \code{2}.
#' @param algorithm String. Specifies the estimation algorithm to use. Default
#'   is \code{"fullrank"}.
#'
#'   Available options include:
#'   \itemize{
#'     \item \code{"fullrank"}: Variational inference with a full-rank multivariate normal approximation.
#'     \item \code{"meanfield"}: Variational inference with independent normal distributions (fast but approximate).
#'     \item \code{"sampling"}: Full Bayesian inference via Markov chain Monte Carlo (MCMC).
#'     \item \code{"pathfinder"}: Quasi-variational inference using optimization-based approximations.
#'     \item \code{"laplace"}: Laplace approximation around the posterior mode.
#'     \item \code{"fixed_param"}: No inference; generates samples using fixed parameter values (e.g., for simulation).
#'   }
#'
#' @param backend Character string specifying which backend to use for fitting
#'   the Stan model. Options are \code{"rstan"} (the default) or
#'   \code{"cmdstanr"}.
#'
#'   Although \code{"rstan"} is the default, it is recommended to use
#'   \code{"cmdstanr"} for better performance and parallelization. To use
#'   \code{cmdstanr}, install the package from GitHub. For more details, see:
#'   \url{https://mc-stan.org/cmdstanr/}
#' @param seed Numeric. The seed for the random number generation for
#'   reproducibility. Default is \code{NA}.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{\code{models}}{A list of model fit objects used to generate the synthetic data.}
#'   \item{\code{data}}{A data frame containing the synthetic data generated from the fitted models.}
#' }
#'
#' @importFrom stats as.formula
#' 
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
    warmup = floor(iter/2),
    thin = 2,
    algorithm = c("fullrank", "meanfield", "sampling", "pathfinder", "laplace", "fixed_param"),
    backend = c("rstan", "cmdstanr"),
    seed = NA
) {
  

  # MATCH ARGUMENTS ---------------------------------------------------------
  
  noise_method <- match.arg(noise_method)
  predictor_method <- match.arg(predictor_method)
  algorithm <- match.arg(algorithm)
  backend <- match.arg(backend)
  

  # FUNCTION ----------------------------------------------------------------
  
  # If variables is NULL, then get all variables
  if(is.null(variables)) {
    
    variables <- setdiff(colnames(data), covariates)
    
  }
  
  # Create formula functions
  formulas <- setNames(
    lapply(
      variables,
      function(x) {
        
        # If predictor method is set to full
        if (predictor_method == "full") {
          
          predictors <- setdiff(variables, x)
          
        }
        
        # If predictor method is set to random
        if (predictor_method == "random") {
          
          predictors <- sample(setdiff(variables, x), predictor_k)
          
        }
        
        # If predictor method is set to information
        if (predictor_method == "information") {
          
          # Get top normalized MI predictors from other variables
          nmi_values <- normalized_mi_scores(
            x, data, variables
          )
          
          # Sort MI values from largest to smallest
          nmi_sorted <- sort(nmi_values, decreasing = TRUE)
          
          # Get top k
          predictors <- names(nmi_sorted)[1:predictor_k]
          
        }
        
        # If predictor method is set to information random mix
        if (predictor_method == "random_information") {
          
          # Get top normalized MI predictors from other variables
          nmi_values <- normalized_mi_scores(
            x, data, variables
          )
          
          # Sort MI values from largest to smallest
          nmi_sorted <- sort(nmi_values, decreasing = TRUE)
          
          # Get top k*5 MI values
          nmi_top <- names(nmi_sorted)[1:(min(predictor_k*5, length(variables)-1))]
          
          predictors <- sample(nmi_top, size = predictor_k, replace = FALSE)
          
        }
        
        # Combine with covariates (if any)
        rhs <- if (!is.null(covariates)) {
          paste0(c(predictors, covariates), collapse = " + ")
        } else {
          paste0(predictors, collapse = " + ")
        }
        
        # Create formula
        as.formula(
          paste0(x, " ~ ", rhs)
        )
        
      }
    ),
    variables)
  
  
  # Randomly sample with replacement ---
  
  # Set arguments
  args <- list(
    data, 
    replace = replace,
    n = n_aug
  )
  
  # Sample data
  sampled_data <- do.call(
    dplyr::slice_sample,
    args
  )
  
  # Inject noise to sampled data
  noisy_data <- add_noise(
    sampled_data,
    variables = variables,
    method = noise_method,
    c = c, 
    p = p, 
    f = f,
    round_int = round_int,
    r = r,
    bounds = noise_bounds
  )
  
  # Get variable families
  variable_families <- assign_families(
    data = data,
    variables = variables,
    manual_families = families
  )
  
  # Create a new progress bar for generating brm models
  pb <- progress::progress_bar$new(
    format = "  Generating Data [:bar] :percent ETA: :eta",
    total  = length(variables),
    clear  = FALSE,
    width  = 100
  )
  
  # Beta shrinkage
  b_shrinkage <- estimate_shrinkage_prior(
    n_obs   = nrow(data),
    n_syn   = n_aug,
    base_sd = base_sd,
    k       = k
  )
  
  # Get predicted values
  brm_models <- setNames(
    lapply(
      variables,
      function(v) {
        
        # Intercept prior
        intercept_prior <- estimate_intercept_prior(
          y        = data[[v]],
          sd_scale = sd_scale
        )
        
        # Fit
        model <- suppressMessages(
          quiet_brm(
            formula   = formulas[[v]],
            data      = data,
            family    = variable_families[[v]],
            prior     = c(b_shrinkage, intercept_prior),
            chains    = chains,
            cores     = if (cores=="auto") {parallel::detectCores()-1} else {cores},
            iter      = iter,
            warmup    = warmup,
            thin      = thin,
            algorithm = algorithm,
            seed      = seed,
            backend   = backend,
            
            # Silence all output messages
            refresh   = 0,
            silent    = 2
          )
        )
        
        # Build newdata
        preds <- noisy_data |> 
          dplyr::select(
            dplyr::all_of(
              all.vars(formulas[[v]])[-1]
            )
          )
        
        # Draw
        pp      <- brms::posterior_predict(model, newdata = preds)
        nrow_pp <- nrow(pp)
        ncol_pp <- ncol(pp)
        rows    <- sample(nrow_pp, ncol_pp, replace = TRUE)
        
        vec     <- as.integer(pp[cbind(rows, seq_len(ncol_pp))])
        
        # Set progress tick
        pb$tick()
        
        # Return list
        return_list <- list(
          model = model,
          value = vec
        )
        
        return(return_list)
        
      }
    ),
    variables
  )
  
  # Extract synthetic predicted data
  synthetic_data <- purrr::set_names(
    purrr::map_dfc(
      brm_models,
      "value"
    ),
    variables
  )
  
  # Extract synthetic models
  synth_models <- purrr::map(
    brm_models,
    "model"
  )
  
  # If there are any covariates, append
  if (!is.null(covariates)) {
    
    cov_data <- noisy_data[, covariates, drop = FALSE]
    
    synthetic_data <- cbind(
      synthetic_data,
      cov_data
    )
    
  }
  
  # Return list
  return_list <- list(
    models = synth_models,
    data = synthetic_data
  )
  
  return(return_list)
  
}