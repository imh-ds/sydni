#' Add Noise to Selected Variables in a Data Frame
#'
#' @description Applies controlled noise to one or more variables in a data
#' frame. It wraps around \code{inject_noise()} and supports multiple methods
#' for injecting noise into numeric, integer, or ordered factor columns.
#'
#' @param data A data frame or tibble.
#' @param variables A tidyselect specification indicating which columns in
#'   \code{data} to modify. If \code{NULL} (default), noise is applied to all
#'   columns.
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
#'   \code{x}.
#'
#' @return A modified data frame of the same structure as \code{data}, with
#'   noise added to the specified variables.
#'
#' @examples
#' df <- data.frame(a = 1:10, b = rnorm(10, 50, 5))
#'
#' # Apply noise to all columns
#' add_noise(df, method = "gaussian", c = 0.05)
#'
#' # Apply discrete noise only to column a
#' add_noise(df, variables = dplyr::all_of("a"), method = "discrete", c = 0.1)
#'
#' @export
add_noise <- function(
    data = .,
    variables = NULL,
    method = c("truncnormal", "gaussian", "discrete", "random"),
    c = "auto",
    round_int = TRUE,
    p = 0.5,
    f = 0.5,
    r = 1,
    bounds = NULL
) {
  
  data |> 
    dplyr::mutate(
      dplyr::across(
        .cols = if (is.null(variables)) {dplyr::everything()} else {variables},
        .fns = function(x) {
          inject_noise(
            x,
            method = method,
            c = c,
            round_int = round_int,
            p = p,
            f = f,
            r = 1,
            bounds = bounds
          )
        }
      )
    )
  
}