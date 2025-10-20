#' Estimate a Prior for the Intercept
#'
#' This function generates a weakly informative normal prior for the intercept
#' in a `brms` model, based on the observed distribution of the outcome
#' variable. The prior is centered at the sample mean and scaled by a multiple
#' of the sample standard deviation to allow flexibility while still applying
#' regularization.
#'
#' @param y Numeric vector. The outcome variable used to estimate the prior.
#' @param sd_scale Numeric. Multiplier for the standard deviation of `y` to set
#'   the prior's spread. A value of 2 approximates a 95% coverage around the
#'   sample mean. Default is 2.
#'
#' @importFrom stats sd
#' 
#' @return A `brms::prior()` object specifying a normal prior for class
#'   `"Intercept"`.
#'
#' @keywords internal
estimate_intercept_prior <- function(
    y,
    sd_scale = 2
) {
  
  # Compute the empirical mean of the outcome
  mean_y <- mean(y, na.rm = TRUE)
  
  # Compute the standard deviation of the outcome
  sd_y <- sd(y, na.rm = TRUE)
  
  # Create and return the brms prior object for the intercept
  do.call(
    brms::prior,
    list(
      
      # Insert computed mean and scaled SD
      prior = bquote(normal(.(mean_y), .(sd_y * sd_scale))),
      
      # Apply to intercept only
      class = "Intercept"
      
    )
  )
  
}
