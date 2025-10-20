#' Quietly Run brms::brm() Without Printing Output
#'
#' This wrapper function runs \code{brms::brm()} while suppressing all printed
#' output, including messages and standard console output (e.g., compilation
#' messages, Stan progress, etc.). It is useful for batch modeling or
#' programmatic workflows where console clutter should be avoided.
#'
#' @param ... Arguments passed directly to \code{brms::brm()}.
#'
#' @return A fitted \code{brmsfit} model object, identical to what
#'   \code{brms::brm()} would return.
#'
#' @importFrom utils capture.output
#'
#' @examples
#' \dontrun{
#' fit <- quiet_brm(
#'   formula = y ~ x,
#'   data = my_data,
#'   family = gaussian()
#' )
#' }
#'
#' @keywords internal
quiet_brm <- function(...) {
  
  # Suppress both messages (warnings, progress) and standard output (compilation info)
  invisible(
    capture.output(
      capture.output(
        fit <- brms::brm(...),  # Fit the model
        type = "message"        # First capture messages (stderr)
      ),
      type = "output"           # Then capture standard output (stdout)
    )
  )
  
  # Return the fitted model object
  fit
  
}
