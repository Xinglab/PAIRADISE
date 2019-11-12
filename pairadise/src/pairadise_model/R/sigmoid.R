#' sigmoid
#'
#'
#' Takes in a vector and applies the sigmoid function elementwise to that vector
#'
#' @name sigmoid
#' @param x : numeric vector
#' @return sigmoid(x)
#' @export 

sigmoid <- function(x){
    1.0 / (1.0 + exp(-x))
}