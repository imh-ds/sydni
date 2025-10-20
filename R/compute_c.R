#' Compute noise scale \code{c} to achieve a target perturbation probability
#'
#' This helper function computes the noise scale \code{c} for the
#' \code{"discrete"} method in \code{inject_noise()}, such that approximately a
#' fraction \code{p} of values are expected to change by at least 1 unit after
#' rounding.
#'
#' @param x A numeric or integer vector. The standard deviation of \code{x} is
#'   used to scale the noise.
#' @param p A numeric value in (0, 1). The target probability that any given
#'   value changes by at least 1 after noise injection and rounding.
#'
#' @return A numeric value \code{c}, which can be passed to
#'   \code{inject_noise()} as the scale for discrete noise.
#'
#' @keywords internal
compute_c_for_p <- function(x, p) {
  
  # x: numeric or integer vector
  # p: desired fraction that should change by >= 1 after rounding
  sigma <- sd(x)
  if (p <= 0 || p >= 1) stop("p must be between 0 and 1, exclusive")
  c_val <- 0.5 / (sigma * (1 - p))
  c_val
  
}

#' Compute noise scale \code{c} such that noise variance equals \code{f *
#' Var(x)}
#'
#' This helper function calculates the appropriate scale \code{c} such that the
#' variance of the injected noise is approximately \code{f} times the variance
#' of the input \code{x}. This is used for the \code{"gaussian"} and
#' \code{"truncnormal"} methods in \code{inject_noise()}.
#'
#' @param f A numeric value in the range (0, 1). Represents the target ratio of
#'   noise variance to the variance of \code{x}, i.e., \code{Var(noise) /
#'   Var(x)}.
#'
#' @return A numeric value \code{c}, suitable for scaling the standard deviation
#'   in noise injection such that \code{Var(noise) = f * Var(x)}.
#'
#' @importFrom stats sd
#'
#' @keywords internal
compute_c_for_f <- function(f) {
  
  # f: desired Var(noise) / Var(x)
  if (f < 0 || f > 1) stop("f should be between 0 and 1")
  sqrt(f)
  
}