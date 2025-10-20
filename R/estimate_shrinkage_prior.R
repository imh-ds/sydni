#' Estimate a Shrinkage Prior for Regression Coefficients
#'
#' This function computes a dynamic prior for regression coefficients (`class =
#' "b"`) in a `brms` model based on the ratio of synthetic to observed data. It
#' returns a `brms::prior()` object with a normal distribution centered at 0 and
#' a standard deviation that shrinks as the proportion of synthetic data
#' increases.
#'
#' @param n_obs Integer. The number of observed (real) data points.
#' @param n_syn Integer. The number of synthetic data points to be added.
#' @param base_sd Numeric. The base standard deviation for the prior when no
#'   synthetic data is used. Default is 10.
#' @param k Numeric. A scaling exponent controlling the rate of shrinkage.
#'   Higher values shrink the prior more aggressively as synthetic data
#'   increases. Default is 1.
#'
#' @return A `brms::prior()` object specifying a normal prior for class `"b"`
#'   coefficients.
#'
#' @keywords internal
estimate_shrinkage_prior <- function(
    n_obs,
    n_syn,
    base_sd = 10,
    k = 1
) {
  
  # Compute the synthetic-to-observed ratio
  r <- n_syn / n_obs
  
  # Adjust the prior standard deviation based on synthetic data ratio
  shrink_sd <- base_sd / (1 + r)^k
  
  # Construct and return the brms prior object for regression coefficients
  do.call(brms::prior, list(
    
    # Insert the evaluated shrink_sd into the prior
    prior = bquote(normal(0, .(shrink_sd))),
    
    # Apply to class "b" (slope coefficients)
    class = "b"
    
  ))
}
