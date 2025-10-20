#' Shrink a Group-Level Correlation Matrix Toward a Global Matrix
#'
#' Applies shrinkage to a group-level correlation matrix by blending it with a
#' global correlation matrix. This is useful for stabilizing correlation
#' estimates in small subgroups by borrowing strength from the overall
#' population.
#'
#' @param cor_group A numeric correlation matrix computed from a subgroup (e.g.,
#'   from upsampled, noise-injected data).
#' @param cor_global A numeric correlation matrix computed from the full dataset
#'   (used as the shrinkage target).
#' @param n_group An integer representing the number of observed rows in the
#'   subgroup.
#' @param shrink_const A positive numeric shrinkage constant (default is 30).
#'   Higher values increase the influence of the global matrix for small-n
#'   groups.
#'
#' @return A numeric correlation matrix of the same dimensions as
#'   \code{cor_group}, representing the shrinkage-weighted average between the
#'   group and global matrices.
#'
#' @examples
#' \dontrun{
#'   cor_shrunk <- shrink_cormat(cor_group, cor_global, n_group = 10)
#' }
#'
#' @keywords internal
shrink_cormat <- function(
    cor_group, 
    cor_global, 
    n_group, 
    shrink_const = 30
) {
  
  lambda <- n_group / (n_group + shrink_const)
  lambda * cor_group + (1 - lambda) * cor_global
  
}
