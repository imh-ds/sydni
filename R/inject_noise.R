#' Inject controlled noise into a numeric, integer, or ordered factor variable
#'
#' Adds a small amount of noise to a single variable using one of several
#' methods: `"discrete"`, `"truncnormal"`, `"gaussian"`, or `"random"`. The
#' function automatically detects whether the input is integer-like, ordered,
#' or continuous and preserves the original range.
#'
#' @param x A vector of type integer, numeric, factor, or ordered factor.
#' @param method Character. One of `"discrete"`, `"truncnormal"`, `"gaussian"`,
#'   or `"random"`.
#' @param c A numeric value in `(0, 1)`. Used to scale the standard deviation of
#'   the noise. Defaults to `"auto"`, in which case it is computed based on
#'   either `p` (for discrete) or `f` (for continuous) inputs.
#' @param round_int Logical. If `TRUE` (default), integer-like inputs are
#'   rounded to the nearest whole number after noise is added.
#' @param p Numeric. Probability that a discrete integer-like value receives a
#'   non-zero perturbation when `c = "auto"`.
#' @param f Numeric. Target ratio of noise variance to data variance for
#'   continuous noise when `c = "auto"`.
#' @param r Numeric. Maximum magnitude of noise when `method = "random"`.
#' @param bounds Optional numeric vector of length 2. If specified, clamps
#'   output to these lower and upper bounds. Otherwise, clamps to the range of
#'   `x`.
#'
#' @details
#' For integer-like or ordered variables, `method = "discrete"` adds noise from
#' `-1`, `0`, or `+1` with probabilities `c / 2`, `1 - c`, and `c / 2`,
#' respectively. For continuous variables, the same method uses uniform noise in
#' `[-c * sd(x), c * sd(x)]`.
#'
#' @return A vector of the same length and compatible type as `x`, with noise
#'   added.
#'
#' @importFrom stats sd rnorm runif
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
  method <- match.arg(method)
  n <- length(x)

  if (c == "auto") {
    c <- if (method == "discrete") {
      compute_c_for_p(x, p)
    } else {
      compute_c_for_f(f)
    }
  }

  is_int <- FALSE

  if (is.ordered(x)) {
    is_int <- TRUE
    x_num <- as.integer(x)
  } else if (is.factor(x)) {
    is_int <- TRUE
    x_num <- as.integer(x)
  } else if (is.integer(x)) {
    is_int <- TRUE
    x_num <- x
  } else if (is.numeric(x)) {
    tol <- sqrt(.Machine$double.eps)
    if (all(abs(stats::na.omit(x) - round(stats::na.omit(x))) < tol)) {
      is_int <- TRUE
      x_num <- as.integer(round(x))
    } else {
      x_num <- x
    }
  } else {
    stop("`x` must be integer, numeric, factor, or ordered factor.")
  }

  sd_x <- stats::sd(as.numeric(x_num))
  if (is.na(sd_x) || sd_x == 0) {
    sd_x <- 1
  }

  eps <- switch(
    method,
    discrete = {
      if (is_int) {
        sample(
          x = c(-1, 0, 1),
          size = n,
          replace = TRUE,
          prob = c(c / 2, 1 - c, c / 2)
        )
      } else {
        stats::runif(n, -c * sd_x, c * sd_x)
      }
    },
    truncnormal = {
      raw <- truncnorm::rtruncnorm(
        n,
        a = -c * sd_x,
        b = c * sd_x,
        mean = 0,
        sd = c * sd_x
      )
      if (is_int) round(raw) else raw
    },
    gaussian = {
      raw <- stats::rnorm(n, mean = 0, sd = c * sd_x)
      if (is_int) round(raw) else raw
    },
    random = {
      raw <- stats::runif(n, min = -r, max = r)
      if (is_int) round(raw) else raw
    }
  )

  x_new <- x_num + eps

  if (is.null(bounds)) {
    lo <- min(x_num, na.rm = TRUE)
    hi <- max(x_num, na.rm = TRUE)
  } else {
    lo <- bounds[1]
    hi <- bounds[2]
  }

  x_new <- pmin(pmax(x_new, lo), hi)
  if (is_int && round_int) {
    x_new <- round(x_new)
  }

  if (is.ordered(x)) {
    level_index <- pmin(pmax(as.integer(x_new), 1L), nlevels(x))
    return(factor(
      levels(x)[level_index],
      levels = levels(x),
      ordered = TRUE
    ))
  }

  if (is.factor(x)) {
    level_index <- pmin(pmax(as.integer(x_new), 1L), nlevels(x))
    return(factor(
      levels(x)[level_index],
      levels = levels(x)
    ))
  }

  if (is.integer(x)) {
    return(as.integer(x_new))
  }

  if (is.numeric(x) && is_int && round_int) {
    return(as.integer(x_new))
  }

  x_new
}
