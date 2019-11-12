#' logit
#'
#'
#' Takes in a vector and applies the logit function elementwise to that vector
#' 
#' @name logit
#' @param x : numeric vector, whose entries should be strictly between 0 and 1
#' @return logit(x)
#' @export 

logit <- function(x){
  if((sum(x >= 1) != 0) | (sum(x <= 0) != 0))
    print('ERROR IN LOGIT: DIVIDING BY ZERO')
  else
    log(x/(1-x))
}