#' Guess an Appropriate Statistical Family for a Response Variable
#'
#' This function inspects a response variable \code{y} and returns a recommended
#' \code{brmsfamily} object for use with \pkg{brms}, based on the variable's
#' type, distribution, and value characteristics.
#'
#' It distinguishes between binary, ordinal, categorical, count, continuous, and
#' proportion outcomes and selects an appropriate family such as
#' \code{bernoulli()}, \code{cumulative()}, \code{categorical()},
#' \code{poisson()}, \code{negbinomial()}, \code{beta()}, or \code{gaussian()}.
#'
#' @param y A vector representing the response variable. Can be a factor,
#'   ordered factor, or numeric vector.
#'
#' @return A \code{brmsfamily} object corresponding to a distribution family,
#'   such as \code{bernoulli("logit")}, \code{gaussian()}, or \code{poisson()}.
#'
#' @importFrom stats gaussian poisson
#' @importFrom brms negbinomial
#'
#' @keywords internal
guess_family <- function(y) {
  
  if (is.factor(y)) {
    n_levels <- nlevels(y)
    
    if (n_levels == 2) return(brms::bernoulli("logit"))
    if (is.ordered(y)) return(brms::cumulative("logit"))
    return(brms::categorical("logit"))
    
  } else if (is.numeric(y)) {
    unique_vals <- length(unique(y))
    is_integer <- all(abs(y - round(y)) < .Machine$double.eps^0.5, na.rm = TRUE)
    
    if (is_integer) {
      max_val <- max(y, na.rm = TRUE)
      min_val <- min(y, na.rm = TRUE)
      
      # Treat small-range integers as ordinal (Likert)
      if (unique_vals >= 3 && unique_vals <= 7 && min_val == 1) {
        return(brms::cumulative("logit"))
      }
      
      # Check for true count-like behavior
      if (all(y >= 0)) {
        if (var(y, na.rm = TRUE) > mean(y, na.rm = TRUE) * 1.5) {
          return(negbinomial())
        } else {
          return(poisson())
        }
      }
    }
    
    if (all(y > 0 & y < 1, na.rm = TRUE)) return(beta())
    
    return(gaussian())
    
  } else {
    
    stop("Unrecognized or unsupported variable type")
    
  }
  
}


#' Assign Statistical Families to Variables
#'
#' Determines appropriate \code{brms} families for a set of variables in a data
#' frame. If a variable is listed in \code{manual_families}, the specified
#' family is used; otherwise, a suitable family is guessed using
#' \code{guess_family()}.
#'
#' @param data A data frame or tibble containing the variables to be modeled.
#' @param variables A character vector of variable names to assign families to.
#'   If \code{NULL} (default), all columns in \code{data} are considered.
#' @param manual_families A named list of families to use manually, where names
#'   correspond to column names in \code{data} and values are \code{brmsfamily}
#'   objects.
#'
#' @return A named list of \code{brmsfamily} objects, one for each specified
#'   variable.
#'
#' @importFrom stats setNames
#'
#' @noRd
#' @keywords internal
assign_families <- function(
    data, 
    variables = NULL, 
    manual_families = list()
) {
  
  if (is.null(variables)) {
    variables <- names(data)
  }
  
  family_list <- setNames(
    lapply(variables, function(var) {
      if (var %in% names(manual_families)) {
        return(manual_families[[var]])  # use manual
      } else {
        return(guess_family(data[[var]]))  # guess it
      }
    }),
    variables)
  
  return(family_list)
}
