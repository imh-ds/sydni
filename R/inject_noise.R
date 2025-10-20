#' Inject controlled noise into a numeric, integer, or ordered factor variable
#'
#' Adds a small amount of noise to a single variable using one of several
#' methods: "discrete", "truncnormal", "gaussian", or "random". The function
#' automatically detects whether the input is an integer, ordered factor, or
#' continuous numeric, and preserves the original value range.
#'
#' @param x A vector of type integer, numeric, or ordered factor.
#' @param method Character. One of "discrete", "truncnormal", "gaussian", or
#'   "random".
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
#' @param bounds Optional numeric vector of length 2. If specified, clamps
#'   output to these lower and upper bounds. Otherwise, clamps to the range of
#'   `x`.
#'
#' @details \strong{Detection of integer-ness:}
#' \enumerate{
#'   \item If `x` is of class integer or an ordered factor, it is treated as integer.
#'   \item If `x` is numeric and all values are effectively whole numbers (within machine tolerance),
#'         it is treated as integer.
#'   \item Otherwise, it is treated as continuous.
#' }
#'
#' After generating noise, the function adds it to `x`, clamps the result to the
#' original or specified range, and rounds if requested. The result preserves
#' the type of `x`.
#'
#' @return A vector of the same length and type as `x`, with noise added.
#' 
#' @importFrom stats sd rnorm runif
#' 
#' @examples
#' x_int  <- sample(1:5, 100, TRUE)
#' inject_noise(x_int, method = "discrete", c = 0.1)
#'
#' x_real <- rnorm(100, 50, 10)
#' inject_noise(x_real, method = "gaussian", c = 0.05, round_int = FALSE)
#'
#' @importFrom truncnorm rtruncnorm
#' @export
inject_noise <- function(
    x,
    method = c("truncnormal", "gaussian", "discrete", "random"),
    c = "auto",
    round_int = TRUE,
    p = 0.5,
    f = 0.5,
    r = 1,
    bounds = NULL
) {
  
  # Set method
  method <- match.arg(method)
  n <- length(x)
  
  # If `c` is set to automate
  if (c == "auto") {
    
    c <- if (method == "discrete") {
      
      compute_c_for_p(x, p)
      
    } else {
      
      compute_c_for_f(f)
      
    }
    
  }
  
  # Detect integer‐ness or ordered factor
  is_int <- FALSE
  
  if (is.ordered(x)) {
    
    is_int <- TRUE
    x_num  <- as.integer(x)
    
  } else if (is.integer(x)) {
    
    is_int <- TRUE
    x_num  <- x
    
  } else if (is.numeric(x)) {
    
    # check if all values are whole numbers (within precision)
    tol <- sqrt(.Machine$double.eps)
    
    if (all(abs(x - round(x)) < tol)) {
      
      is_int <- TRUE
      x_num  <- as.integer(round(x))
      
    } else {
      
      x_num  <- x
      
    }
    
  } else {
    
    stop("`x` must be integer, numeric, or ordered factor.")
    
  }
  
  # Compute base SD
  sd_x <- sd(as.numeric(x_num))
  
  # Generate noise ε
  eps <- switch(
    method,
    
    discrete = {
      # draw uniform noise in [-c·sd, +c·sd]
      raw <- runif(n, -c * sd_x, c * sd_x)
      if (is_int) round(raw) else raw
    },
    
    truncnormal = {
      a <- -c * sd_x
      b <- c * sd_x
      raw <- truncnorm::rtruncnorm(n, a=a, b=b, mean=0, sd=c * sd_x)
      if (is_int) round(raw) else raw
    },
    
    gaussian = {
      raw <- rnorm(n, mean=0, sd=c * sd_x)
      if (is_int) round(raw) else raw
    },
    
    random = {
      raw <- runif(n, min = -r, max = r)
      if (is_int) round(raw) else raw
    },
    
    {  # fallback
      warning("`method` not recognized; defaulting to 'truncnormal'")
      a <- -c * sd_x
      b <- c * sd_x
      raw <- truncnorm::rtruncnorm(n, a=a, b=b, mean=0, sd=c * sd_x)
      if (is_int) round(raw) else raw
    }
    
  )
  
  # Apply, clamp, and (optionally) round
  x_new <- x_num + eps
  
  if (is.null(bounds)) {
    
    lo    <- min(x_num)
    hi    <- max(x_num)
    
  } else {
    
    lo <- bounds[1]
    hi <- bounds[2]
    
  }

  x_new <- pmin(pmax(x_new, lo), hi)
  if (is_int && round_int) x_new <- round(x_new)
  
  # Return in original type
  if (is.ordered(x)) {
    
    levels <- sort(unique(x_num))
    factor(x_new,
           levels = levels,
           ordered = TRUE)
    
  } else if (is.integer(x)) {
    
    as.integer(x_new)
    
  } else if (is.numeric(x) && is_int && round_int) {
    
    as.integer(x_new)
    
  } else {
    
    x_new
    
  }
}