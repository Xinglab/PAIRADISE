#!/usr/bin/env Rscript
########################################
## PAIRADISE power simulation
## Author: Levon Demirdjian
## Last Updated: 10/16/2019 
## Description: This script implements
## a PAIRADISE simulation based on the 
## model parameters estimated on the 
## Geuvadis CEU population, with the 
## number of reads manually specified
########################################

library(PAIRADISE)

## Define sigmoid function
sigmoid <- function(x){
  1/(1 + exp(-x))
}

## Define logit function
logit <- function(x){
  log(x/(1 - x))
}

## Read in model parameters
filename   <- 'Results_CEU/CEU_parameters_filtered.txt'
parameters <- read.table(file = filename, header = TRUE)

## Save individual parameters
sigma1_dist <- parameters$sigma1
sigma2_dist <- parameters$sigma2
sigma_dist  <- parameters$sigma
delta_dist  <- parameters$delta
mu_dist     <- parameters$mu
lI_dist     <- parameters$lI
lS_dist     <- parameters$lS
nEvents_CEU <- nrow(parameters)

## Zero out the middle 50% of the empirical estimates of delta to generate null and alternative cases
delta_adj <- delta_dist
delta_adj[delta_dist >= quantile(delta_dist, 0.25) & delta_dist <= quantile(delta_dist, 0.75)] <- 0

## Compute 25th, 50th, 75th quantiles of sigma1, sigma2, sigma
sigma1 <- quantile(x = sigma1_dist, probs = c(0.25, 0.50, 0.75))
sigma2 <- quantile(x = sigma2_dist, probs = c(0.25, 0.50, 0.75))
sigma  <- quantile(x = sigma_dist,  probs = c(0.25, 0.50, 0.75))

## Read in user inputs (nReplicates: positive integer, var_in: low, medium, high)
args        <- commandArgs(trailingOnly = TRUE)
K           <- as.numeric(args[1]) ## K = number of replicates
var_in      <- as.numeric(args[2]) ## 1 = low, 2 = medium, 3 = high
n_reads     <- as.numeric(args[3]) ## n_reads = number of RNA-seq reads
sigma_index <- var_in

## Randomly shuffle indices (with replacement) for later sampling
shuffled_ind <- sample(x = nEvents_CEU, size = 50000, replace = TRUE)

## Generate table of simulated data
nEvents  <- 5000
sim_data <- matrix(0, nrow = nEvents, ncol = 12)

## Name the columns
colnames(sim_data) <- c('ID', 'I1', 'S1', 'I2', 'S2', 'lI', 'lS', 'delta', 'mu', 'logit_psi1', 'logit_psi2', 'alpha')
k <- 1
for(i in 1:nEvents){
  
  ## Randomly sample a set of parameters (jointly)
  delta_ind <- delta_adj[shuffled_ind[i]]
  mu_ind    <- mu_dist[shuffled_ind[i]]
  lI_ind    <- lI_dist[shuffled_ind[i]]
  lS_ind    <- lS_dist[shuffled_ind[i]]
  
  ## Set variance
  sigma1_ind <- sigma1[sigma_index]
  sigma2_ind <- sigma2[sigma_index]
  sigma_ind  <- sigma[sigma_index]
  
  ## Generate logit(psi1), logit(psi2), alpha
  alpha      <- rnorm(n = K, mean = mu_ind, sd = sigma_ind)
  logit_psi1 <- rnorm(n = K, mean = alpha,  sd = sigma1_ind)
  logit_psi2 <- rnorm(n = K, mean = alpha + delta_ind, sd = sigma2_ind)
  
  ## Length-normalize psi1 and psi2
  psi1_adj <- lI_ind * sigmoid(logit_psi1) / (lI_ind * sigmoid(logit_psi1) + lS_ind * (1 - sigmoid(logit_psi1)))
  psi2_adj <- lI_ind * sigmoid(logit_psi2) / (lI_ind * sigmoid(logit_psi2) + lS_ind * (1 - sigmoid(logit_psi2)))

  ## Generate I1, I2, S1, S2
  I1 <- rbinom(n = K, size = n_reads, prob = psi1_adj)
  I2 <- rbinom(n = K, size = n_reads, prob = psi2_adj)
  S1 <- n_reads - I1
  S2 <- n_reads - I2
  
  ## Save data to matrix
  sim_data[k, 1] <- paste('Event_', i, sep = '')
  sim_data[k, 2] <- paste(I1, collapse = ',')
  sim_data[k, 3] <- paste(S1, collapse = ',')
  sim_data[k, 4] <- paste(I2, collapse = ',')
  sim_data[k, 5] <- paste(S2, collapse = ',')
  sim_data[k, 6] <- lI_ind
  sim_data[k, 7] <- lS_ind
  
  ## Save parameters to matrix
  sim_data[k, 8]  <- delta_ind
  sim_data[k, 9]  <- mu_ind
  sim_data[k, 10] <- paste(logit_psi1, collapse = ',')
  sim_data[k, 11] <- paste(logit_psi2, collapse = ',')
  sim_data[k, 12] <- paste(alpha, collapse = ',')
  k <- k + 1

  if(k > nEvents)
    break
  
}

## Write data to table
sim_data_filename <- paste('Simulated_Data_CEU/SimData_', K, '_replicates_', var_in, '_variance_', n_reads, '_reads.txt', sep = '')
write.table(x = sim_data, file = sim_data_filename, row.names = FALSE, col.names = TRUE, quote = FALSE)

## Run PAIRADISE on the simulated data
pairadise_results <- pairadise(data.frame(sim_data[,1:7]))

## Save workspace
save.image(paste('Results_CEU/SimData_', K, '_replicates_', var_in, '_variance_', n_reads, '_reads.RData', sep = ''))
