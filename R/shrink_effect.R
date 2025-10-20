#' Shrink Effect Sizes Toward Marginal Distributions
#'
#' This function applies shrinkage to a dataset by partially replacing observed values
#' with draws from the empirical marginal distributions of each column. The effect is to
#' reduce the overall magnitude of associations or effects, mimicking a form of noise
#' injection while preserving univariate structure.
#'
#' @param data A data frame or tibble containing only numeric or ordinal values. Each
#'   column is treated as a discrete variable (e.g., Likert scale).
#' @param c A numeric value in (0, 1). Controls the strength of the shrinkage.
#'   Smaller values produce stronger shrinkage. Internally, the shrinkage intensity
#'   is computed as \code{lambda = 1 - sqrt(c)}.
#'
#' @return A data frame of the same shape as \code{data}, with values partially
#'   resampled from their marginal distributions and rounded to the nearest integer.
#'   All values are clamped to the range 1–5.
#' 
#' @importFrom stats rnorm sd
#'
#' @keywords internal
shrink_effect <- function(
  data,
  c
) {
  
  n <- nrow(data)
  p <- ncol(data)
  
  # 1) compute lambda
  lambda <- 1 - sqrt(c)
  
  # 2) empirical marginals for each item
  marginals <- lapply(data, function(x) {
    tab <- table(x)/length(x)
    probs <- as.numeric(tab); names(probs) <- names(tab)
    probs
  })
  
  # 3) sample noise U
  U <- sapply(1:p, function(j) {
    sample(
      x    = as.numeric(names(marginals[[j]])),
      size = n,
      replace = TRUE,
      prob = marginals[[j]]
    )
  })
  U <- as.data.frame(U)
  colnames(U) <- colnames(data)
  
  # 4) mix & round
  data_shrunk <- as.data.frame(
    round((1 - lambda) * as.matrix(data) + lambda * as.matrix(U))
  )
  # enforce 1–5 bounds (in case of jitter to 0 or 6)
  data_shrunk[] <- lapply(data_shrunk, function(x) pmin(pmax(x,1),5))
  
  return(data_shrunk)
  
}


#' Shrink Effect Sizes with Sign-Preserving Noise
#'
#' Applies a shrinkage transformation to the input data while preserving the
#' general directionality (sign structure) of effects. Noise is added with mean
#' zero and matched scale, allowing signal degradation without flipping overall
#' relationships.
#'
#' Unlike \code{shrink_effect()}, this version adds Gaussian noise rather than
#' resampling from marginal distributions.
#'
#' @param data A data frame or matrix containing only numeric values. Typically
#'   represents ordinal or continuous Likert-style data.
#' @param c A numeric value in (0, 1). Controls shrinkage strength. Higher
#'   values retain more of the original signal; smaller values yield stronger
#'   shrinkage. Internally converted to \code{lambda = 1 - sqrt(c)}.
#' @param max_lambda A numeric upper limit for \code{lambda}. Prevents
#'   over-shrinking. Default is \code{0.6}.
#'
#' @return A data frame of the same dimensions as \code{data}, with Gaussian
#'   noise applied and values clamped to the range [1, 5] after rounding.
#'
#' @noRd
#' @keywords internal
shrink_effect_sign_preserving <- function(
    data, 
    c, 
    max_lambda = 0.6
) {
  
  n <- nrow(data)
  p <- ncol(data)
  lambda <- min(1 - sqrt(c), max_lambda)
  
  # generate noise: mean-zero and matched scale
  noise <- matrix(rnorm(n * p, mean = 0, sd = sd(as.matrix(data))), nrow = n)
  
  # apply shrinkage
  data_shrunk <- as.data.frame((1 - lambda) * as.matrix(data) + lambda * noise)
  
  # round and bound
  data_shrunk[] <- lapply(data_shrunk, function(x) round(pmin(pmax(x, 1), 5)))
  
  return(data_shrunk)
  
}
