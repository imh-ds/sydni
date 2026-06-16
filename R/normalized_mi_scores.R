#' Compute Normalized Mutual Information Scores
#'
#' This function computes the normalized mutual information (MI) between a target variable
#' and a set of other variables. Normalization is done by dividing the raw MI by the maximum
#' entropy of the two variables, resulting in a score between 0 and 1.
#'
#' @param target A character string indicating the name of the target variable in \code{data}.
#' @param data A data frame containing all relevant variables. Variables should be discrete
#'   (e.g., factors or integers).
#' @param variables A character vector specifying which columns in \code{data} to include
#'   in the comparison. The target should be a subset of these.
#'
#' @return A named numeric vector of normalized MI scores, one for each variable in
#'   \code{variables} excluding the \code{target}. Each score lies in the range [0, 1],
#'   with higher values indicating stronger dependency with the target.
#'
#' @keywords internal
normalized_mi_scores <- function(
    target,
    data, 
    variables
) {
  discretize_if_needed <- function(x) {
    if (is.factor(x) || is.ordered(x)) {
      return(as.integer(x))
    }

    if (is.numeric(x) || is.integer(x)) {
      ux <- unique(stats::na.omit(x))
      if (length(ux) <= 10) {
        return(as.integer(factor(x, exclude = NULL)))
      }

      breaks <- unique(stats::quantile(x, probs = seq(0, 1, length.out = 6), na.rm = TRUE))
      if (length(breaks) < 3) {
        return(as.integer(factor(x, exclude = NULL)))
      }
      return(as.integer(cut(x, breaks = breaks, include.lowest = TRUE, ordered_result = TRUE)))
    }

    stop("Mutual-information predictor selection supports numeric and factor columns only.")
  }

  others <- setdiff(variables, target)
  target_x <- discretize_if_needed(data[[target]])
  
  sapply(others, function(var) {
    other_x <- discretize_if_needed(data[[var]])
    mi <- infotheo::mutinformation(target_x, other_x)
    h1 <- infotheo::entropy(target_x)
    h2 <- infotheo::entropy(other_x)
    mi / max(h1, h2)  # or sqrt(h1 * h2)
  })
  
}
