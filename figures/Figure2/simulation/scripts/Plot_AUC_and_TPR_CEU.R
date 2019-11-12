#!/usr/bin/env Rscript

##########################################
## PAIRADISE simulation AUC and TPR plots
## Author: Levon Demirdjian
## Last Updated: 08/19/2019 
## Description: This script plots the 
## AUC and TPR curves for the simulation 
## results of PAIRADISE
##########################################

## Input:
## 1) table_loc: Directory where the simulation tables are located

## Directory where the simulation tables are located
args      <- commandArgs(TRUE)
table_loc <- args[1]

## Read in the data
AUC_filenames <- paste(table_loc, 'AUC_', c('Low', 'Medium', 'High'), '_Variance_CEU.txt', sep = '')
TPR_filenames <- paste(table_loc, 'TPR_', c('Low', 'Medium', 'High'), '_Variance_CEU.txt', sep = '')

AUC_list <- list()
TPR_list <- list()
for(i in 1:3){
  AUC_list[[i]] <- read.table(file = AUC_filenames[i], header = TRUE)
  TPR_list[[i]] <- read.table(file = TPR_filenames[i], header = TRUE)
}

######################################
## Make AUC plots for each variance ##
######################################

pdf('AUC_Curves_CEU.pdf', width = 10, height = 5)
par(mfrow = c(1, 3))
par(pty = 's')

##################
## Low Variance ##
##################

plot(1:5, AUC_list[[1]][1,], col = 'black', ylim = c(0.5, 1), xaxt = 'n',
     main = 'Low Variance', xlab = 'Number of Replicates', ylab = 'Area Under the Curve (AUC)')
lines(1:5, AUC_list[[1]][1,], type = 'l', col = 'black')
axis(1, at = 1:5, labels = c(3, 5, 10 ,20 ,50))

points(1:5, AUC_list[[1]][2,], col = 'blue')
lines(1:5, AUC_list[[1]][2,], type = 'l', col = 'blue', lty = 2)

points(1:5, AUC_list[[1]][3,], col = 'purple')
lines(1:5, AUC_list[[1]][3,], type = 'l', col = 'purple', lty = 2)

points(1:5, AUC_list[[1]][4,], col = 'green')
lines(1:5, AUC_list[[1]][4,], type = 'l', col = 'green', lty = 2)

points(1:5, AUC_list[[1]][5,], col = 'red')
lines(1:5, AUC_list[[1]][5,], type = 'l', col = 'red', lty = 2)


legend(3, 0.9, legend = c("PAIRADISE", "Paired T", "Paired Wilcoxon", "rMATS Paired", "Fisher's Test"),
       col = c("black", "blue", "purple", "red", "green"), lty = c(1, 2, 2, 2, 2), cex = 0.8)

#####################
## Medium Variance ##
#####################

plot(1:5, AUC_list[[2]][1,], col = 'black', ylim = c(0.5, 1), xaxt = 'n',
     main = 'Medium Variance', xlab = 'Number of Replicates', ylab = 'Area Under the Curve (AUC)')
lines(1:5, AUC_list[[2]][1,], type = 'l', col = 'black')
axis(1, at = 1:5, labels = c(3, 5, 10 ,20 ,50))

points(1:5, AUC_list[[2]][2,], col = 'blue')
lines(1:5, AUC_list[[2]][2,], type = 'l', col = 'blue', lty = 2)

points(1:5, AUC_list[[2]][3,], col = 'purple')
lines(1:5, AUC_list[[2]][3,], type = 'l', col = 'purple', lty = 2)

points(1:5, AUC_list[[2]][4,], col = 'green')
lines(1:5, AUC_list[[2]][4,], type = 'l', col = 'green', lty = 2)

points(1:5, AUC_list[[2]][5,], col = 'red')
lines(1:5, AUC_list[[2]][5,], type = 'l', col = 'red', lty = 2)


###################
## High Variance ##
###################

plot(1:5, AUC_list[[3]][1,], col = 'black', ylim = c(0.5, 1), xaxt = 'n',
     main = 'High Variance', xlab = 'Number of Replicates', ylab = 'Area Under the Curve (AUC)')
lines(1:5, AUC_list[[3]][1,], type = 'l', col = 'black')
axis(1, at = 1:5, labels = c(3, 5, 10 ,20 ,50))

points(1:5, AUC_list[[3]][2,], col = 'blue')
lines(1:5, AUC_list[[3]][2,], type = 'l', col = 'blue', lty = 2)

points(1:5, AUC_list[[3]][3,], col = 'purple')
lines(1:5, AUC_list[[3]][3,], type = 'l', col = 'purple', lty = 2)

points(1:5, AUC_list[[3]][4,], col = 'green')
lines(1:5, AUC_list[[3]][4,], type = 'l', col = 'green', lty = 2)

points(1:5, AUC_list[[3]][5,], col = 'red')
lines(1:5, AUC_list[[3]][5,], type = 'l', col = 'red', lty = 2)

dev.off()

######################################
## Make TPR plots for each variance ##
######################################

pdf('TPR_Curves_CEU.pdf', width = 10, height = 5)
par(mfrow = c(1, 3))
par(pty = 's')

##################
## Low Variance ##
##################

plot(1:5, TPR_list[[1]][1,], col = 'black', ylim = c(0, 1), xaxt = 'n',
     main = 'Low Variance', xlab = 'Number of Replicates', ylab = 'TPR at 5% FPR')
lines(1:5, TPR_list[[1]][1,], type = 'l', col = 'black')
axis(1, at = 1:5, labels = c(3, 5, 10 ,20 ,50))

points(1:5, TPR_list[[1]][2,], col = 'blue')
lines(1:5, TPR_list[[1]][2,], type = 'l', col = 'blue', lty = 2)

points(1:5, TPR_list[[1]][3,], col = 'purple')
lines(1:5, TPR_list[[1]][3,], type = 'l', col = 'purple', lty = 2)

points(1:5, TPR_list[[1]][4,], col = 'green')
lines(1:5, TPR_list[[1]][4,], type = 'l', col = 'green', lty = 2)

points(1:5, TPR_list[[1]][5,], col = 'red')
lines(1:5, TPR_list[[1]][5,], type = 'l', col = 'red', lty = 2)

#####################
## Medium Variance ##
#####################

plot(1:5, TPR_list[[2]][1,], col = 'black', ylim = c(0, 1), xaxt = 'n',
     main = 'Medium Variance', xlab = 'Number of Replicates', ylab = 'TPR at 5% FPR')
lines(1:5, TPR_list[[2]][1,], type = 'l', col = 'black')
axis(1, at = 1:5, labels = c(3, 5, 10 ,20 ,50))

points(1:5, TPR_list[[2]][2,], col = 'blue')
lines(1:5, TPR_list[[2]][2,], type = 'l', col = 'blue', lty = 2)

points(1:5, TPR_list[[2]][3,], col = 'purple')
lines(1:5, TPR_list[[2]][3,], type = 'l', col = 'purple', lty = 2)

points(1:5, TPR_list[[2]][4,], col = 'green')
lines(1:5, TPR_list[[2]][4,], type = 'l', col = 'green', lty = 2)

points(1:5, TPR_list[[2]][5,], col = 'red')
lines(1:5, TPR_list[[2]][5,], type = 'l', col = 'red', lty = 2)


###################
## High Variance ##
###################

plot(1:5, TPR_list[[3]][1,], col = 'black', ylim = c(0, 1), xaxt = 'n',
     main = 'High Variance', xlab = 'Number of Replicates', ylab = 'TPR at 5% FPR')
lines(1:5, TPR_list[[3]][1,], type = 'l', col = 'black')
axis(1, at = 1:5, labels = c(3, 5, 10 ,20 ,50))

points(1:5, TPR_list[[3]][2,], col = 'blue')
lines(1:5, TPR_list[[3]][2,], type = 'l', col = 'blue', lty = 2)

points(1:5, TPR_list[[3]][3,], col = 'purple')
lines(1:5, TPR_list[[3]][3,], type = 'l', col = 'purple', lty = 2)

points(1:5, TPR_list[[3]][4,], col = 'green')
lines(1:5, TPR_list[[3]][4,], type = 'l', col = 'green', lty = 2)

points(1:5, TPR_list[[3]][5,], col = 'red')
lines(1:5, TPR_list[[3]][5,], type = 'l', col = 'red', lty = 2)

dev.off()

