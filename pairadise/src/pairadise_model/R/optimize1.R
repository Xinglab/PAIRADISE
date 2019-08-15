#' optimize1
#'
#'
#' Used internally in PAIRADISE to compute the MLEs of delta, mu, sigma1, sigma2, sigma
#'
#' @name optimize1
#' @param x Numeric vector such that x = (sigma1, sigma2, sigma, mu, delta) if equal.variance = FALSE, and x = (sigma1, sigma, mu, delta) if equal.variance = TRUE. sigma1, sigma2, sigma must be positive
#' @param M Number of replicates for the current exon.
#' @param I1 Exon inclusion counts for group 1. Positive integers.
#' @param S1 Exon skipping counts for group 1. Positive integers.
#' @param I2 Exon inclusion counts for group 2. Positive integers.
#' @param S2 Exon skipping counts for group 2. Positive integers.
#' @param l.iI Effective length of inclusion isoform. Positive integer.
#' @param l.iS Effective length of skipping isoform. Positive integer.
#' @param logit.psi1 Numeric vector with values of logit psi1.
#' @param logit.psi2 Numeric vector with values of logit psi2.
#' @param alpha Numeric vector with values of alpha.
#' @param equal.variance Are the group variances assumed equal? Default value is FALSE.

#' @export


optimize1 <- function(x, M, I1, S1, I2, S2, l.iI, l.iS, logit.psi1, logit.psi2, alpha, equal.variance = FALSE){

  numerator1 <- l.iI * l.iS * sigmoid(logit.psi1) * (sigmoid(logit.psi1) - 1) * (I1 + S1)
  numerator2 <- l.iI * l.iS * sigmoid(logit.psi2) * (sigmoid(logit.psi2) - 1) * (I2 + S2)

  denominator1 <- (l.iI * sigmoid(logit.psi1) + l.iS * (1 - sigmoid(logit.psi1)))^2
  denominator2 <- (l.iI * sigmoid(logit.psi2) + l.iS * (1 - sigmoid(logit.psi2)))^2

  if(equal.variance == TRUE){
  
    a1 <- M * (2 * log(x[1]) + log(x[2]))
    a2 <- sum((logit.psi1 - alpha)^2) / (2*(x[1]^2))
    a3 <- sum((logit.psi2 - alpha - x[4])^2) / (2*(x[1]^2))
    a4 <- sum((alpha - x[3])^2) / (2*(x[2]^2))
  
    b1 <- (1/x[1]^4) * ((1/(x[1]^2)) - numerator2 / denominator2)
    b2 <- (1/(x[1]^2)) - (numerator1 / denominator1)
    b3 <- (1/x[1]^4) + ((2/(x[1]^2)) + (1/(x[2]^2))) * ((numerator2 / denominator2) - (1/(x[1]^2)))

  }else{
    
    a1 <- M * (log(x[1]) + log(x[2]) + log(x[3]))
    a2 <- sum((logit.psi1 - alpha)^2) / (2*(x[1]^2))
    a3 <- sum((logit.psi2 - alpha - x[5])^2) / (2*(x[2]^2))
    a4 <- sum((alpha - x[4])^2) / (2*(x[3]^2))
    
    b1 <- (1/x[1]^4) * ((1/(x[2]^2)) - numerator2 / denominator2)
    b2 <- (1/(x[1]^2)) - (numerator1 / denominator1)
    b3 <- (1/x[2]^4) + ((1/(x[1]^2)) + (1/(x[2]^2)) + (1/(x[3]^2))) * ((numerator2 / denominator2) - (1/(x[2]^2)))
    
  }
    
  result <- a1 + a2 + a3 + a4 + 0.25 * sum(log((b1 + b2 * b3)^2))
  result

}
