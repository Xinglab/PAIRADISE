#' optimize2
#'
#'
#' Used internally in PAIRADISE to compute the MLEs of logit(psi1), logit(psi2), alpha 
#'
#' @name optimize2
#' @param x Numeric vector such that x = (logit(psi1), logit(psi2), alpha)
#' @param k Index representing current replicate number.
#' @param I1 Exon inclusion counts for group 1. Positive integers. 
#' @param S1 Exon skipping counts for group 1. Positive integers.
#' @param I2 Exon inclusion counts for group 2. Positive integers.
#' @param S2 Exon skipping counts for group 2. Positive integers.
#' @param l.iI Effective length of inclusion isoform. Positive integer.
#' @param l.iS Effective length of skipping isoform. Positive integer.
#' @param mu Parameter mu.
#' @param delta Parameter delta.
#' @param s1 Group 1 standard deviation. Positive number.
#' @param s2 Group 2 standard deviation. Positive number.
#' @param s Overall standard deviation. Positive number.
#' 
#' @export 


optimize2 <- function(x, k, I1, S1, I2, S2, l.iI, l.iS, delta, mu, s1, s2, s){
          
  a1 <- ((x[1] - x[3])^2) / (2*(s1^2))
  a2 <- ((x[2] - x[3] - delta)^2) / (2*(s2^2))
  a3 <- ((x[3] - mu)^2) / (2*(s^2))
          
  b1 <- I1[k] * log((l.iI * sigmoid(x[1])) / (l.iI * sigmoid(x[1]) + l.iS * (1 - sigmoid(x[1]))))
  b2 <- S1[k] * log((l.iS * (1 - sigmoid(x[1]))) / (l.iI * sigmoid(x[1]) + l.iS * (1 - sigmoid(x[1]))))
  b3 <- I2[k] * log((l.iI * sigmoid(x[2])) / (l.iI * sigmoid(x[2]) + l.iS * (1 - sigmoid(x[2]))))
  b4 <- S2[k] * log((l.iS * (1 - sigmoid(x[2]))) / (l.iI * sigmoid(x[2]) + l.iS * (1 - sigmoid(x[2]))))

  result <- a1 + a2 + a3 - (b1 + b2 + b3 + b4)
  result
        
}