
#' Add Backticks to Variable Names
#'
#' This function adds backticks around variable names to ensure proper formatting in formulas.
#'
#' @param x A vector of strings representing variable names.
#'
#' @return A vector of strings where each element is enclosed in backticks.
#' @export
#'
#' @examples
#' add_backticks(c("age", "weight", "height"))
#' # Output: "`age`", "`weight`", "`height`"
add_backticks = function(x) {
  paste0("`", x, "`")
}


#' Construct Formula for Independent Variables
#'
#' This function creates a formula string for the predictor variables, with each variable enclosed in backticks and concatenated with "+".
#'
#' @param x A vector of strings representing predictor variable names.
#'
#' @return A single string that concatenates variable names with "+" for use in a model formula.
#' @export
#'
#' @examples
#' x_lm_formula(c("age", "weight", "height"))
#' # Output: "`age` + `weight` + `height`"
x_lm_formula = function(x) {
  paste(add_backticks(x), collapse = "+")
}


#' Build Linear Model Formula
#'
#' This function creates a complete linear model formula using the provided response variable and predictor variables.
#'
#' @param x A vector of predictor variable names.
#' @param y A single response variable name.
#'
#' @return A formula object suitable for use in linear models such as \code{lm()}.
#' @export
#'
#' @examples
#' build_lm_formula(c("age", "weight", "height"), "outcome")
#' # Output: outcome ~ `age` + `weight` + `height`
build_lm_formula = function(x, y){
  if (length(y)>1){
    stop("y needs to be just one variable")
  }
  as.formula(
    paste0("`",y,"`", " ~ ", x_lm_formula(x))
  )
}
