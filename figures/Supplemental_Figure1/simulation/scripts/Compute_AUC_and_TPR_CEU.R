#!/usr/bin/env Rscript
##########################################
## PAIRADISE simulation AUC and TPR plots
## Author: Levon Demirdjian
## Last Updated: 10/16/2019 
## Description: This script computes the 
## AUC and TPR curves for the power 
## simulation results of PAIRADISE 
##########################################

## Input:
## 1) res_loc: Directory where the simulation is being performed

library(AUC)      
library(PAIRADISE)

## Working directory where the simulation is being performed
args    <- commandArgs(TRUE)
res_loc <- args[1]

## Compute AUC and TPR curves
variance_vector  <- c(1, 2, 3)
replicate_vector <- c(3, 5, 10, 20, 50)
reads_vector     <- c(10, 20, 50, 100, 1000)

AUC_list         <- lapply(1:3, function(x) matrix(0, 5, 5))
names(AUC_list)  <- c('Low_Variance', 'Medium_Variance', 'High_Variance')

TPR_list         <- lapply(1:3, function(x) matrix(0, 5, 5))
names(TPR_list)  <- c('Low_Variance', 'Medium_Variance', 'High_Variance')

## Label rows and columns for the different methods
for(i in 1:3){
  rownames(AUC_list[[i]]) <- reads_vector
  rownames(TPR_list[[i]]) <- reads_vector
  colnames(AUC_list[[i]]) <- replicate_vector
  colnames(TPR_list[[i]]) <- replicate_vector
}

for(read_ind in 1:5){
  for(var_ind in 1:3){
    for(rep_ind in 1:5){
    
      ## Read in simulation results
      cur_var <- variance_vector[var_ind]
      K       <- replicate_vector[rep_ind]
      load(paste(res_loc, 'Results_CEU/SimData_', K,  '_replicates_', cur_var, '_variance_', reads_vector[read_ind], '_reads.RData', sep = ''))
   
      ## Read in p-values and hypotheses
      pairadise_pvals <- as.numeric(pairadise_results$raw.pvalues)
      true_deltas     <- as.numeric(sim_data[, 8])
      true_classes    <- 1 - 1 * (true_deltas == 0)
    
      ## Compute PAIRADISE AUC
      roc_pairadise <- roc(1 - pairadise_pvals, factor(true_classes))
      auc_pairadise <- auc(roc_pairadise)
    
      ## Compute PAIRADISE TPR
      upper_ind <- which(roc_pairadise$fpr > 0.05)[1]
      lower_ind <- upper_ind - 1
      fpr_lower <- roc_pairadise$fpr[lower_ind]
      fpr_upper <- roc_pairadise$fpr[upper_ind]
    
      ## Linear interpolation
      y <- (0.05 - fpr_lower) / (fpr_upper - fpr_lower)
      tpr_pairadise <- roc_pairadise$tpr[lower_ind] + y * (roc_pairadise$tpr[upper_ind] - roc_pairadise$tpr[lower_ind])
    
      ## Save results into AUC and TPR list
      AUC_list[[var_ind]][read_ind, rep_ind] <- auc_pairadise
      TPR_list[[var_ind]][read_ind, rep_ind] <- tpr_pairadise
   
    } 
  }
}

## Save results into three tables (one for each variance setting)
for(i in 1:3){
  write.table(x = AUC_list[[i]], file = paste(res_loc, '/Results_CEU/',  'AUC_', names(AUC_list)[i], '_CEU.txt', sep = ''), quote = FALSE, col.names = NA)
}

## Save results into three tables (one for each variance setting)
for(i in 1:3){
  write.table(x = TPR_list[[i]], file = paste(res_loc, '/Results_CEU/',  'TPR_', names(TPR_list)[i], '_CEU.txt', sep = ''), quote = FALSE, col.names = NA)
}

