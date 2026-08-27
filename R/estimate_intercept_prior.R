#' Estimate a Prior for the Intercept
#'
#' This function generates a weakly informative normal prior for the intercept
#' in a Gaussian `brms` model, based on the observed distribution of the
#' outcome variable. The prior is centered at the sample mean and scaled by a
#' multiple of the sample standard deviation. Non-Gaussian outcomes return
#' `NULL` so that `brms` uses its family-appropriate default intercept prior.
#'
#' @param y Outcome vector used to estimate a Gaussian-model prior.
#' @param sd_scale Numeric. Multiplier for the standard deviation of `y` to set
#'   the prior's spread. A value of 2 approximates a 95% coverage around the
#'   sample mean. Default is 2.
#' @param family Optional model family. When supplied and not Gaussian, no
#'   empirical normal intercept prior is created.
#'
#' @importFrom stats sd
#' 
#' @return A `brms::prior()` object specifying a normal prior for class
#'   `"Intercept"`.
#'
#' @keywords internal
estimate_intercept_prior <- function(
    y,
    sd_scale = 2,
    family = NULL
) {
  if (
    is.factor(y) ||
      !is.numeric(y) ||
      (!is.null(family) && !identical(family$family, "gaussian"))
  ) {
    return(NULL)
  }
  
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
