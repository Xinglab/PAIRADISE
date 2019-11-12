#!/usr/bin/env Rscript
##########################################
## PAIRADISE simulation AUC and TPR plots
## Author: Levon Demirdjian
## Last Updated: 08/19/2019 
## Description: This script computes the 
## AUC and TPR curves for the simulation 
## results of PAIRADISE
##########################################

## Input:
## 1) res_loc: Directory where the simulation is being performed

library(AUC)      
library(PAIRADISE)

## Working directory where the simulation is being performed
args    <- commandArgs(TRUE)
res_loc <- args[1]

## Function to compute p-values using the other (competing) methods
alternative_tests <- function(sim_data){
  
  ## Save read counts and effective lengths
  I1 <- sim_data[,2]
  S1 <- sim_data[,3]
  I2 <- sim_data[,4]
  S2 <- sim_data[,5]
  lI <- as.numeric(sim_data[,6])
  lS <- as.numeric(sim_data[,7])
  n  <- nrow(sim_data)
  
  paired_t_pvals <- rep(0, n)
  wilcox_pvals   <- rep(0, n)
  fisher_pvals   <- rep(0, n)
  
  for(i in 1:n){
    
    ## Unroll the read counts
    I1_ind <- as.numeric(strsplit(x = I1[i], split = ',')[[1]])
    S1_ind <- as.numeric(strsplit(x = S1[i], split = ',')[[1]])
    I2_ind <- as.numeric(strsplit(x = I2[i], split = ',')[[1]])
    S2_ind <- as.numeric(strsplit(x = S2[i], split = ',')[[1]])
    
    ## Compute naive estimates of psi
    psi1_naive <- lS[i] * I1_ind / (lS[i] * I1_ind + lI[i] * S1_ind)
    psi2_naive <- lS[i] * I2_ind / (lS[i] * I2_ind + lI[i] * S2_ind)
    
    if(is.na(psi1_naive) | is.na(psi2_naive))
      next

    ## Note: the t-test and Wilcoxon t-tests might not work for some events 
    ## since the counts can be constant. If they don't work, we assign a 
    ## p-value of 1.
    
    ## Compute t-test p-values
    t_try <- try(t.test(x = psi1_naive, y = psi2_naive, paired = TRUE, alternative = 'two.sided'))
    if("try-error" %in% class(t_try)){
      paired_t_pvals[i] <- 1
    }else{
      paired_t_pvals[i] <- t.test(x = psi1_naive, y = psi2_naive, paired = TRUE, alternative = 'two.sided')$p.value
    }
    #t_results <- t.test(x = psi1_naive, y = psi2_naive, paired = TRUE, alternative = 'two.sided')
    
    
    ## Compute Wilcoxon p-values
    w_try <- try(wilcox.test(x = psi1_naive, y = psi2_naive, paired = TRUE, alternative = 'two.sided'))
    if("try-error" %in% class(w_try)){
      wilcox_pvals[i] <- 1
    }else{
      wilcox_pvals[i] <- wilcox.test(x = psi1_naive, y = psi2_naive, paired = TRUE, alternative = 'two.sided')$p.value
    }
    #wilcox_results  <- wilcox.test(x = psi1_naive, y = psi2_naive, paired = TRUE, alternative = 'two.sided')

    ## Compute Fisher p-values
    fisher_pvals_all <- rep(0, K)
    for(k in 1:K){
      fisher_table        <- rbind(c(I1_ind[k], S1_ind[k]), c(I2_ind[k], S2_ind[k]))
      fisher_results      <- fisher.test(x = fisher_table, alternative = 'two.sided')
      fisher_pvals_all[k] <- fisher_results$p.value
    }
    
    ## Combine Fisher p-values
    fisher_stat     <- -2 * sum(log(fisher_pvals_all))
    fisher_pvals[i] <- dchisq(x = fisher_stat, df = 2 * K)
    
  }
  
  ## Return a list with all of the p-values
  list(paired_t_pvals = paired_t_pvals, wilcox_pvals = wilcox_pvals, fisher_pvals = fisher_pvals)
  
}

## Compute AUC and TPR curves
variance_vector  <- c(1, 2, 3)
replicate_vector <- c(3, 5, 10, 20, 50)

AUC_list         <- lapply(1:3, function(x) matrix(0, 5, 5))
names(AUC_list)  <- c('Low_Variance', 'Medium_Variance', 'High_Variance')

TPR_list         <- lapply(1:3, function(x) matrix(0, 5, 5))
names(TPR_list)  <- c('Low_Variance', 'Medium_Variance', 'High_Variance')

## Label rows and columns for the different methods
for(i in 1:3){
  rownames(AUC_list[[i]]) <- c('PAIRADISE', 't-test', 'Wilcox', 'Fisher', 'rMATS')
  rownames(TPR_list[[i]]) <- c('PAIRADISE', 't-test', 'Wilcox', 'Fisher', 'rMATS')
  colnames(AUC_list[[i]]) <- c(3, 5, 10 ,20 ,50)
  colnames(TPR_list[[i]]) <- c(3, 5, 10 ,20 ,50)
}

for(var_ind in 1:3){
  for(rep_ind in 1:5){
    
    ## Read in simulation results
    cur_var <- variance_vector[var_ind]
    K       <- replicate_vector[rep_ind]
    load(paste(res_loc, 'Results_CEU/SimData_', K,  '_replicates_', cur_var, '_variance.RData', sep = ''))
   
    ## Read in rMATS results
    rMATS_floc    <- paste(res_loc, 'Results_CEU/rMATS/', sep = '')
    rMATS_fname   <- paste(rMATS_floc, K, '_replicates_', cur_var, '_variance/rMATS_Result_P_', K, '_replicates_', cur_var, '_variance.txt', sep = '')
    rMATS_results <- read.table(rMATS_fname, header = TRUE, sep = '\t')
 
    ## Read in p-values and hypotheses
    pairadise_pvals <- as.numeric(pairadise_results$raw.pvalues)
    true_deltas     <- as.numeric(sim_data[, 8])
    true_classes    <- 1 - 1 * (true_deltas == 0)

    ## Compute p-values for all of the other methods
    alt_pvals <- alternative_tests(sim_data)
    
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
    
    ## Compute t-test AUC
    roc_t_test <- roc(1 - alt_pvals$paired_t_pvals, factor(true_classes))
    auc_t_test <- auc(roc_t_test)
    
    ## Compute t-test TPR
    upper_ind <- which(roc_t_test$fpr > 0.05)[1]
    lower_ind <- upper_ind - 1
    fpr_lower <- roc_t_test$fpr[lower_ind]
    fpr_upper <- roc_t_test$fpr[upper_ind]
    
    y <- (0.05 - fpr_lower) / (fpr_upper - fpr_lower)
    tpr_t_test <- roc_t_test$tpr[lower_ind] + y * (roc_t_test$tpr[upper_ind] - roc_t_test$tpr[lower_ind])
    
    ## Compute Wilcoxon AUC
    roc_wilcox_test <- roc(1 - alt_pvals$wilcox_pvals, factor(true_classes))
    auc_wilcox_test <- auc(roc_wilcox_test)
    
    ## Compute Wilcoxon TPR
    upper_ind <- which(roc_wilcox_test$fpr > 0.05)[1]
    lower_ind <- upper_ind - 1
    fpr_lower <- roc_wilcox_test$fpr[lower_ind]
    fpr_upper <- roc_wilcox_test$fpr[upper_ind]
    
    y <- (0.05 - fpr_lower) / (fpr_upper - fpr_lower)
    tpr_wilcox_test <- roc_wilcox_test$tpr[lower_ind] + y * (roc_wilcox_test$tpr[upper_ind] - roc_wilcox_test$tpr[lower_ind])
    
    ## Compute Fisher AUC
    roc_fisher_test <- roc(1 - alt_pvals$fisher_pvals, factor(true_classes))
    auc_fisher_test <- auc(roc_fisher_test)
    
    ## Compute Fisher TPR
    upper_ind <- which(roc_fisher_test$fpr > 0.05)[1]
    lower_ind <- upper_ind - 1
    fpr_lower <- roc_fisher_test$fpr[lower_ind]
    fpr_upper <- roc_fisher_test$fpr[upper_ind]
    
    y <- (0.05 - fpr_lower) / (fpr_upper - fpr_lower)
    tpr_fisher_test <- roc_fisher_test$tpr[lower_ind] + y * (roc_fisher_test$tpr[upper_ind] - roc_fisher_test$tpr[lower_ind])

    ## Compute rMATS-paired AUC
    roc_rMATS <- roc(1 - rMATS_results$PValue, factor(true_classes))
    auc_rMATS <- auc(roc_rMATS)
    
    ## Compute rMATS-paired TPR
    upper_ind <- which(roc_rMATS$fpr > 0.05)[1]
    lower_ind <- upper_ind - 1
    fpr_lower <- roc_rMATS$fpr[lower_ind]
    fpr_upper <- roc_rMATS$fpr[upper_ind]
    
    y <- (0.05 - fpr_lower) / (fpr_upper - fpr_lower)
    tpr_rMATS <- roc_rMATS$tpr[lower_ind] + y * (roc_rMATS$tpr[upper_ind] - roc_rMATS$tpr[lower_ind])
    
    ## Save results into AUC and TPR list
    AUC_list[[var_ind]][,rep_ind] <- c(auc_pairadise, auc_t_test, auc_wilcox_test, auc_fisher_test, auc_rMATS)
    TPR_list[[var_ind]][,rep_ind] <- c(tpr_pairadise, tpr_t_test, tpr_wilcox_test, tpr_fisher_test, tpr_rMATS)
    
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

