##########################################
## PAIRADISE analysis of CEU population
## Author: Levon Demirdjian
## Last Updated: 08/16/2019 
## Description: This script implements
## PAIRADISE on the CEU population then 
## saves the resulting model parameters
##########################################

library(PAIRADISE)

## Define input file
data_loc <- 'data/CEU_SE_allexons_count.txt'
my.data  <- read.table(data_loc, header = TRUE)[,1:7]
colnames(my.data) <- c('Exon_ID', 'I1', 'S1', 'I2', 'S2', 'lI', 'lS')

## Filter out events where both mean(psi1_hat) and mean(psi2_hat) are > 0.95 or < 0.05
n   <- nrow(my.data)
tag <- rep(0, n)
inv <- rep(0, n)
for(i in 1:n){
  
  ## Save read counts and effective lengths
  I1 <- as.numeric(strsplit(as.character(my.data[i, 2]), split = ',')[[1]])
  S1 <- as.numeric(strsplit(as.character(my.data[i, 3]), split = ',')[[1]])
  I2 <- as.numeric(strsplit(as.character(my.data[i, 4]), split = ',')[[1]])
  S2 <- as.numeric(strsplit(as.character(my.data[i, 5]), split = ',')[[1]])
  lI <- my.data[i,6]
  lS <- my.data[i,7]
  
  valid <- (I1 + S1 != 0) & (I2 + S2 != 0)
  if(is.na(sum(valid))){
    inv[i] <- 1
    next
  }
  
  if(sum(valid) == 0){
    inv[i] <- 1
    next
  }
  
  ## Compute naive estimates of psi
  psi1_naive <- lS * I1[valid] / (lS * I1[valid] + lI * S1[valid])
  psi2_naive <- lS * I2[valid] / (lS * I2[valid] + lI * S2[valid])
  
  if(any(is.na(psi1_naive)) | any(is.na(psi2_naive))){
    inv[i] <- 1
    next
  }
  
  ## Mark if event falls outside of our filtering psi range
  psi1_good <- (mean(psi1_naive) >= 0.05) & (mean(psi1_naive) <= 0.95)
  psi2_good <- (mean(psi2_naive) >= 0.05) & (mean(psi2_naive) <= 0.95)

  if(!(psi1_good | psi2_good))
    tag[i] <- 1
    
}

## Focus on valid subset of the data
my.data_refined <- my.data[tag == 0,]

## Run PAIRADISE
results <- pairadise(my.data[,1:7])

## Save workspace
system('mkdir Results_CEU')
system('mkdir Simulated_Data_CEU')
save.image('Results_CEU/CEU_PAIRADISE_Results_filtered.RData')

## Save parameters for simulations
mu     <- as.numeric(sapply(results$param.unconstrained, function(x) x[1,1]))
sigma1 <- as.numeric(sapply(results$param.unconstrained, function(x) x[2,1]))
sigma2 <- as.numeric(sapply(results$param.unconstrained, function(x) x[3,1]))
sigma  <- as.numeric(sapply(results$param.unconstrained, function(x) x[4,1]))
delta  <- as.numeric(sapply(results$param.unconstrained, function(x) x[5,1]))
lI     <- my.data[,6]
lS     <- my.data[,7]

## Separate the n1 and n2 values with commas
n1 <- sapply(1:nrow(my.data), function(x) paste(as.character(results$I1[[x]] + results$S1[[x]]), collapse = ','))
n2 <- sapply(1:nrow(my.data), function(x) paste(as.character(results$I2[[x]] + results$S2[[x]]), collapse = ','))

## Save data into a table
sim_parameters <- data.frame(sigma1 = sigma1, sigma2 = sigma2, sigma = sigma, delta = delta, mu = mu, lI = lI, lS = lS, n1 = n1, n2 = n2)

## Write table into a .txt file
write.table(x = sim_parameters, file = 'Results_CEU/CEU_parameters_filtered.txt', quote = FALSE, col.names = TRUE, row.names = FALSE)


