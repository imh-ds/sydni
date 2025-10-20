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
  
  others <- setdiff(variables, target)
  
  sapply(others, function(var) {
    mi <- infotheo::mutinformation(data[[target]], data[[var]])
    h1 <- infotheo::entropy(data[[target]])
    h2 <- infotheo::entropy(data[[var]])
    mi / max(h1, h2)  # or sqrt(h1 * h2)
  })
  
}