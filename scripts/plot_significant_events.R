#############################################
## Plot ASAS events identified by PAIRADISE
## Author: Levon Demirdjian
## Last updated: 07/23/2019
## Email: levondem@gmail.com
## Description: This script generates a pdf
## plot showing all of the ASAS events
## identified by PAIRADISE
#############################################

########################
## Script information ##
########################

## Usage: Rscript plot/plot_all_sig_SNPs.R plot_input_data.txt 1 FIN
##
## Inputs:
## 1. plot_input_data.txt (tab-delimited): This is the FDR adjusted (filtered) output of the PAIRADISE statistical package.
##                   ID, reference allele inclusion counts, skipping counts, alternative allele inclusion counts, skipping counts
##                   incLen, skpLen, pval, IncLevel1, IncLevel2, AvgIncLevel1, AvgIncLevel2, IncLevelDifference,
##                   AvgTotalCount1, AvgTotalCount2, FDR
## 2. eff_len_I: 1 or 2: effective length of the exon inclusion isoform
## 3. sample_name: Sample name for the title of the resulting pdf.

## Load required packages
library(binom)
library(ggplot2)

## Read in data
args      <- commandArgs(trailingOnly = TRUE)
my.data   <- read.table(args[1], header = TRUE, stringsAsFactors = FALSE)  ## nSig_events rows
pval      <- formatC(my.data$pval, format = "e", digits = 2)
qval      <- formatC(my.data$qval, format = "e", digits = 2)
eff_len_I <- as.numeric(args[2])

gene_name  <- substr(my.data$ExonID, regexpr('-\"', my.data$ExonID, fixed = TRUE) + 2, regexpr('\"-', my.data$ExonID, fixed = TRUE) - 1)
snp_name   <- substr(my.data$ExonID, regexpr("_rs", my.data$ExonID, ignore.case = FALSE) + 1, regexpr("\\_[^\\_]*$", my.data$ExonID) - 1)
chrom      <- substr(my.data$ExonID, regexpr("chr", my.data$ExonID), regexpr("(", my.data$ExonID, fixed = TRUE) - 1)
pos        <- substr(my.data$ExonID, regexpr(":", my.data$ExonID) + 1, regexpr("_chr", my.data$ExonID, ignore.case = FALSE) - 1)
samples    <- strsplit(my.data$SampleName, split = ',')

allele_loc  <- regexpr("\\_[^\\_]*$", my.data$ExonID)
allele_info <- sapply(1:nrow(my.data), function(x) substr(my.data$ExonID[x], allele_loc[x] + 1, nchar(my.data$ExonID[x]))[[1]])

ref_allele   <- substr(x = allele_info, 1, regexpr('|', allele_info, fixed = TRUE) - 1)
alt_allele   <- substr(x = allele_info, regexpr('|', allele_info, fixed = TRUE) + 1, nchar(allele_info))
allele_freqs <- my.data$AF

## Round AFs to 3 significant digits
#af_1 <- round(as.numeric(substr(x = allele_freqs, 3, regexpr('|', allele_freqs, fixed = TRUE) - 1)), 3)
#af_2 <- round(as.numeric(substr(x = allele_freqs, regexpr('|', allele_freqs, fixed = TRUE) + 3, nchar(allele_freqs))), 3)

#al_end       <- substr(x = allele_freqs, regexpr('|', allele_freqs, fixed = TRUE) + 1, regexpr('|', allele_freqs, fixed = TRUE) + 2)
#allele_freqs <- paste(substr(x = allele_freqs, 1, 2), af_1, '|', al_end, af_2, sep = '')

## Sort everything by p values
ind     <- order(as.numeric(pval))
pval    <- pval[ind]
qval    <- qval[ind]
my.data <- my.data[ind, ]

samples      <- samples[ind]
gene_name    <- gene_name[ind]
snp_name     <- snp_name[ind]
chrom        <- chrom[ind]
pos          <- pos[ind]
ref_allele   <- ref_allele[ind]
alt_allele   <- alt_allele[ind]
allele_freqs <- allele_freqs[ind]

pdf_name <- as.character(args[3])
pdf(paste('plots/', pdf_name, '_significant_events.pdf', sep = ''), width = 8, height = 8)

for(i in 1:nrow(my.data)){

  ## Current sample IDs
  cur_sample <- samples[[i]]
  #cur_sample <- 1:length(as.numeric(strsplit(toString(my.data$IJC_REF[i]),',')[[1]]))

  ## The inputs are allele1 inclusion counts, skipping counts, allele 2 inclusion counts and skipping counts
  I_ref <- as.numeric(strsplit(toString(my.data$IJC_REF[i]),',')[[1]])
  S_ref <- as.numeric(strsplit(toString(my.data$SJC_REF[i]),',')[[1]])
  I_alt <- as.numeric(strsplit(toString(my.data$IJC_ALT[i]),',')[[1]])
  S_alt <- as.numeric(strsplit(toString(my.data$SJC_ALT[i]),',')[[1]])

  ## Remove invalid replicates
  valid1 <- which(I_ref + S_ref != 0)
  valid2 <- which(I_alt + S_alt != 0)
  valid  <- intersect(valid1, valid2)

  I_ref <- I_ref[valid]
  S_ref <- S_ref[valid]
  I_alt <- I_alt[valid]
  S_alt <- S_alt[valid]

  cur_sample <- cur_sample[valid]

  ## Naive psi's
  psi_ref <- I_ref / (I_ref + S_ref)
  psi_alt <- I_alt / (I_alt + S_alt)
  my_ordering <- order(abs(psi_ref - psi_alt))

  ## Reorder the psi's
  I_ref <- I_ref[my_ordering]
  S_ref <- S_ref[my_ordering]
  I_alt <- I_alt[my_ordering]
  S_alt <- S_alt[my_ordering]

  cur_sample <- cur_sample[my_ordering]

  ## Sample ids (we don't want a clogged y-axis)
  if(length(cur_sample) <= 10){
    id <- cur_sample
  }else{
    id <- seq(from = 1, to = length(cur_sample), 1)
  }

  ## Create data frame with reference allele confidence intervals
  df_1        <- cbind(id, binom.confint(I_ref, I_ref + eff_len_I * S_ref, methods = 'lrt')[c(4,5,6)])
  df_1$allele <- ref_allele[i]
  names(df_1) <- c("ID", "rate", "lower", "upper", "Allele")

  ## Create data frame with alternative allele confidence intervals
  df_2        <- cbind(id, binom.confint(I_alt, I_alt + eff_len_I * S_alt, methods = 'lrt')[c(4,5,6)])
  df_2$allele <- alt_allele[i]
  names(df_2) <- c("ID", "rate", "lower", "upper", "Allele")

  ## Combine data frames into one big data frame
  df <- rbind(df_1, df_2)

  df$Allele <- factor(df$Allele, levels = c(ref_allele[i], alt_allele[i]))

  ## Save MAF frequency
  maf_freq <- allele_freqs[i]
  maf_freq <- gsub(maf_freq, pattern = '|', replacement = ', ', fixed = TRUE)
  maf_freq <- gsub(maf_freq, pattern = ':', replacement = ': ', fixed = TRUE)  
# maf_freq <- round(as.numeric(strsplit(alt_freq[i], ',')[[1]]), digits = 4)
  # if(length(maf_freq) > 1){
  #   maf_freq <- paste(maf_freq, collapse = ',')
  # }

  ## Plot the results
  pp <- ggplot(data = df, aes(x = ID, y = rate, ymin = lower, ymax = upper, colour = Allele)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(position = position_dodge(width = 0.5), width = 0.1) +
    scale_colour_manual(values = c("blue", "red")) +
    theme_bw() +
    coord_flip() +
    ylab("Exon Inclusion Level") +
    xlab("Sample") +
    ylim(0, 1) +
    labs(title = paste(snp_name[i], '\n Gene = ', gene_name[i], ' \n', pos[i], ' \n p = ', pval[i], ' \n q = ', qval[i], ' \n AF = ', maf_freq, sep = '')) +
    theme(panel.grid.major.x = element_line(colour = "grey", linetype = "dashed"),
          panel.grid.major.y = element_line(colour = "grey", linetype = "dashed"),
          panel.grid.minor.y = element_line(colour = "grey", linetype = "dashed"),
          plot.title = element_text(hjust = 0.5))

  print(pp)

}

dev.off()

