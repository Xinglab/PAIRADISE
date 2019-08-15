#' loglikelihood
#'
#'
#' Used internally in PAIRADISE to compute the log-likelihood function 
#'
#' @name loglikelihood
#' @param M Number of replicates for the current exon. Positive integer. 
#' @param I1 Exon inclusion counts for group 1. Positive integers. 
#' @param S1 Exon skipping counts for group 1. Positive integers.
#' @param I2 Exon inclusion counts for group 2. Positive integers.
#' @param S2 Exon skipping counts for group 2. Positive integers.
#' @param l.iI Effective length of inclusion isoform. Positive integer.
#' @param l.iS Effective length of skipping isoform. Positive integer.
#' @param logit.psi1 Numeric vector with values of logit psi1.
#' @param logit.psi2 Numeric vector with values of logit psi2.
#' @param alpha Numeric vector with values of alpha.
#' @param s1 Group 1 standard deviation. Positive number.
#' @param s2 Group 2 standard deviation. Positive number.
#' @param s Overall standard deviation. Positive number.
#' @param mu Parameter mu.
#' @param delta Parameter delta.
#' @return log likelihood value at input.
#' 
#' @export 



loglikelihood <- function(M, I1, S1, I2, S2, l.iI, l.iS, logit.psi1, logit.psi2, alpha, s1, s2, s, mu, delta){

  numerator1 <- l.iI * l.iS * sigmoid(logit.psi1) * (sigmoid(logit.psi1) - 1) * (I1 + S1)
  numerator2 <- l.iI * l.iS * sigmoid(logit.psi2) * (sigmoid(logit.psi2) - 1) * (I2 + S2)

  denominator1 <- (l.iI * sigmoid(logit.psi1) + l.iS * (1 - sigmoid(logit.psi1)))^2
  denominator2 <- (l.iI * sigmoid(logit.psi2) + l.iS * (1 - sigmoid(logit.psi2)))^2

  a1 <- M * (log(s1) + log(s2) + log(s))
  a2 <- sum((logit.psi1 - alpha)^2) / (2*(s1^2))
  a3 <- sum((logit.psi2 - alpha - delta)^2) / (2*(s2^2))
  a4 <- sum((alpha - mu)^2) / (2*(s^2))

  b1 <- (1/s1^4) * ((1/(s2^2)) - numerator2 / denominator2)
  b2 <- (1/(s1^2)) - (numerator1 / denominator1)
  b3 <- (1/s2^4) + ((1/(s1^2)) + (1/(s2^2)) + (1/(s^2))) * ((numerator2 / denominator2) - (1/(s2^2)))

  c1 <- sum(I1 * log((l.iI * sigmoid(logit.psi1)) / (l.iI * sigmoid(logit.psi1) + l.iS * (1 - sigmoid(logit.psi1)))))
  c2 <- sum(S1 * log((l.iS * (1 - sigmoid(logit.psi1))) / (l.iI * sigmoid(logit.psi1) + l.iS * (1 - sigmoid(logit.psi1)))))
  c3 <- sum(I2 * log((l.iI * sigmoid(logit.psi2)) / (l.iI * sigmoid(logit.psi2) + l.iS * (1 - sigmoid(logit.psi2)))))
  c4 <- sum(S2 * log((l.iS * (1 - sigmoid(logit.psi2))) / (l.iI * sigmoid(logit.psi2) + l.iS * (1 - sigmoid(logit.psi2)))))

  result <- -(a1 + a2 + a3 + a4 + 0.25 * sum(log((b1 + b2 * b3)^2))) + c1 + c2 + c3 + c4
  result

}

# loglikelihood <- function(M, I1, S1, I2, S2, l.iI, l.iS, logit.psi1, logit.psi2, alpha, s1, s2, s, mu, delta){
# 
#   numerator1 <- l.iI * sigmoid(logit.psi1)
#   numerator2 <- l.iS * (1 - sigmoid(logit.psi1))
#   numerator3 <- l.iI * sigmoid(logit.psi2)
#   numerator4 <- l.iS * (1 - sigmoid(logit.psi2))
# 
#   denominator1 <- (l.iI * sigmoid(logit.psi1) + l.iS * (1 - sigmoid(logit.psi1)))
#   denominator2 <- (l.iI * sigmoid(logit.psi2) + l.iS * (1 - sigmoid(logit.psi2)))
# 
#   a1 <- sum(I1 * log(numerator1 / denominator1))
#   a2 <- sum(S1 * log(numerator2 / denominator1))
#   a3 <- sum(I2 * log(numerator3 / denominator2))
#   a4 <- sum(S2 * log(numerator4 / denominator2))
# 
#   b1 <- M * (log(s1) + log(s2) + log(s))
#   b2 <- sum((logit.psi1 - alpha)^2) / (2*(s1^2))
#   b3 <- sum((logit.psi2 - alpha - delta)^2) / (2*(s2^2))
#   b4 <- sum((alpha - mu)^2) / (2*(s^2))
# 
#   result <- a1 + a2 + a3 + a4 - b1 - b2 - b3 - b4
#   result
# 
# }