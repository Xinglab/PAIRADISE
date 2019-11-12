###############################################
## Run the PAIRADISE statistical model
## Author: Levon Demirdjian
## Last Updated: 07/27/2019
## Description: This script runs the PAIRADISE
## statistical model and creates several output 
## files.
###############################################

## Load the PAIRADISE statistical model
library('PAIRADISE')

###############################################
## Step 1: Run PAIRADISE and get raw results ##
###############################################

## Read in command line inputs
args <- commandArgs(TRUE)
args

## args: 1) Input file ('../output/ASAS.SNP.SE.JunctionReadsOnly.byPair.filtered.txt')
##       2) Result file ('pairadise_result.txt')
##       3) Equal variance? ('FALSE')
##       4) FDR cutoff (0.10)

## Read in the data and run PAIRADISE
my.data <- read.table(args[1], header = TRUE, colClasses = c(rep("character", 5), rep("numeric", 2), rep("NULL", 10)))
results <- pairadise(my.data, numCluster = 8, equal.variance = args[3])

## Write the event IDs along with their raw p-values
raw_name <- gsub(x = args[2], pattern = '.txt', replacement = '_raw.txt')
res.tab  <- cbind(AS_SNP_ID = results$exonID, pvals = results$raw.pvalues)
write.table(data.frame("ID" = rownames(res.tab), res.tab), file = raw_name, row.names = FALSE, quote = FALSE)

## Save sample ID's
sample_ids <- read.table(args[1], header = TRUE, colClasses = c(rep("NULL", 7), rep("character", 1), rep("NULL", 9)))[,1]

## Save alleles and their frequencies
allele_freq <- read.table(args[1], header = TRUE, colClasses = c(rep("NULL", 15), rep("character", 1), rep("NULL", 1)))[,1]
allele_loc  <- regexpr("\\_[^\\_]*$", my.data$ExonID)
allele_info <- sapply(1:nrow(my.data), function(x) substr(my.data$ExonID[x], allele_loc[x] + 1, nchar(my.data$ExonID[x]))[[1]])

## Save R workspace
workspace_name <- gsub(x = args[2], pattern = '.txt', replacement = '.RData')
save.image(workspace_name)

################################################
## Step 2: Append data information to results ##
################################################

input_IDs <- my.data$ExonID
SNP_IDs   <- results$exonID

## Match PAIRADISE result SNP_IDs to the raw input rows
matched_index <- match(input_IDs, SNP_IDs)
pvals_matched <- results$raw.pval[matched_index]

## Write a new table with the p-values appended to the input
#write.table(data.frame(my.data, Samples = sample_ids, AF = allele_freq, pval = pvals_matched), file = args[2], quote = FALSE, row.names = FALSE, sep = '\t')

####################################################
## Step 3: Append inclusion and count information ##
####################################################

## Average total counts
AvgTotalCount1 <- sapply(1:length(SNP_IDs), function(x) mean(results$I1[[x]] + results$S1[[x]]))
AvgTotalCount2 <- sapply(1:length(SNP_IDs), function(x) mean(results$I2[[x]] + results$S2[[x]]))

## Get the effective lengths for these events
lI <- my.data$incLen[matched_index]
lS <- my.data$skpLen[matched_index]

## Inclusion levels
I1 <- my.data$IJC_REF[matched_index]
S1 <- my.data$SJC_REF[matched_index]
I2 <- my.data$IJC_ALT[matched_index]
S2 <- my.data$SJC_ALT[matched_index]

psi1 <- sapply(1:length(SNP_IDs), function(x) paste(round((as.numeric(strsplit(I1[[x]], ',')[[1]])/lI[x]) / ((as.numeric(strsplit(I1[[x]], ',')[[1]])/lI[x]) + (as.numeric(strsplit(S1[[x]], ',')[[1]])/lS[x])), 3), collapse = ','))
psi2 <- sapply(1:length(SNP_IDs), function(x) paste(round((as.numeric(strsplit(I2[[x]], ',')[[1]])/lI[x]) / ((as.numeric(strsplit(I2[[x]], ',')[[1]])/lI[x]) + (as.numeric(strsplit(S2[[x]], ',')[[1]])/lS[x])), 3), collapse = ','))

## Average inclusion levels
valid_ind1 <- sapply(1:length(SNP_IDs), function(x) which(!is.na(as.numeric(strsplit(psi1[x], ',')[[1]]))))
valid_ind2 <- sapply(1:length(SNP_IDs), function(x) which(!is.na(as.numeric(strsplit(psi2[x], ',')[[1]]))))
valid_ind  <- sapply(1:length(SNP_IDs), function(x) intersect(valid_ind1[[x]], valid_ind2[[x]]))

mean_psi1  <- round(sapply(1:length(SNP_IDs), function(x) mean(as.numeric(strsplit(psi1[x], ',')[[1]])[valid_ind[[x]]])), 3)
mean_psi2  <- round(sapply(1:length(SNP_IDs), function(x) mean(as.numeric(strsplit(psi2[x], ',')[[1]])[valid_ind[[x]]])), 3)

## Inclusion level difference
mean_psi_diff <- mean_psi1 - mean_psi2

## FDR adjustment of p-values
qvals <- p.adjust(pvals_matched, method = "BH")

## Write output
final_dataset <- data.frame(my.data, pval = pvals_matched, qval = qvals, IncLevel1 = psi1, IncLevel2 = psi2,
                            AvgIncLevel1 = mean_psi1, AvgIncLevel2 = mean_psi2, IncLevelDifference = mean_psi_diff,
                            AvgTotalCount1 = AvgTotalCount1, AvgTotalCount2 = AvgTotalCount2,
                            SampleName = sample_ids,
                            RefAltAllele = allele_info, AF = allele_freq)

final_name <- gsub(x = args[2], pattern = '.txt', replacement = '_totalcount.txt')
#write.table(final_dataset, file = final_name, quote = FALSE, row.names = FALSE, sep = '\t')

## Save events satisfying the following two criteria: 
## 1) There are >= 10 average total counts in both groups
## 2) At least one group has its average psi value in (0.05, 0.95)

filter_1     <- which(AvgTotalCount1 >= 10 & AvgTotalCount2 >= 10)
filter_2     <- which(((0.05 <= mean_psi1) & (mean_psi1 <= 0.95)) | ((0.05 <= mean_psi2) & (mean_psi2 <= 0.95)))
valid_events <- intersect(filter_1, filter_2)

filtered_name          <- gsub(x = args[2], pattern = '.txt', replacement = '_filtered.txt')
final_dataset_filtered <- final_dataset[valid_events,]

## Recalculate the qvalues using only the filtered events
qvals <- p.adjust(final_dataset_filtered$pval, method = "BH")
final_dataset_filtered$qval <- qvals

write.table(final_dataset_filtered, file = filtered_name, quote = FALSE, row.names = FALSE, sep = '\t')

## Save events significant at FDR <= fdr_thresh
fdr_thresh        <- as.numeric(args[4])
fdr_adjusted_name <- gsub(x = args[2], pattern = '.txt', replacement = paste('_filtered_FDR', fdr_thresh * 100, '.txt', sep = ''))
write.table(final_dataset_filtered[final_dataset_filtered$qval <= fdr_thresh, ], file = fdr_adjusted_name, quote = FALSE, row.names = FALSE, sep = '\t')

