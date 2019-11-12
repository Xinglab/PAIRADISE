#' clean.data
#'
#'
#' Removes missing data and invalid pairs from the matched pair data to be analyzed by PAIRADISE.
#'
#' @name clean.data
#' @param my.data  Data frame containing grouped data to be analyzed. 
#' @return The function clean.data returns a list containing the following entries:
#' \item{I1}{Group 1 isoform 1 counts for each replicate.}
#' \item{S1}{Group 1 isoform 2 counts for each replicate.}
#' \item{I2}{Group 2 isoform 1 counts for each replicate.}
#' \item{S2}{Group 2 isoform 2 counts for each replicate.}
#' \item{length_I}{Effective lengths of isoform 1.}
#' \item{length_S}{Effective lengths of isoform 2.}
#' \item{exonList}{IDs of the exons/events.}
#' \item{nExon}{Number of exons/events.}
#' \item{M}{Vector containing the number of replicates per exon/event.}
#' @details The data frame has 7 columns, arranged as follows:
#' Column 1 contains the ID of the exons/events.
#' Column 2 contains counts of isoform 1 corresponding to the first group. 
#' Column 3 contains counts of isoform 2 corresponding to the first group. 
#' Column 4 contains counts of isoform 1 corresponding to the second group. 
#' Column 5 contains counts of isoform 2 corresponding to the second group. 
#' Replicates in columns 2-5 should be separated by commas, e.g. 1623,432,6 for three replicates.
#' Column 6 contains the effective length of isoform 1.
#' Column 7 contains the effective length of isoform 2.
#' @export 

clean.data <- function(my.data){
    
  ## Unpack data list
  data.list <- load.data(my.data)
  I1.raw <- data.list$I1.raw
  S1.raw <- data.list$S1.raw
  I2.raw <- data.list$I2.raw
  S2.raw <- data.list$S2.raw
  length_I.raw <- data.list$length_I.raw
  length_S.raw <- data.list$length_S.raw
  exonList.raw <- data.list$exonList.raw
  nExon.raw <- data.list$nExon.raw
  M <- data.list$M
  
  ## Clean missing/invalid data
  miss_total <- c()
  missIndex <- 1
  for(iExon in 1:nExon.raw){
    
    ## Find total number of NA values.
    I1.miss <- sum(is.na(I1.raw[[iExon]]))
    S1.miss <- sum(is.na(S1.raw[[iExon]]))
    I2.miss <- sum(is.na(I2.raw[[iExon]]))
    S2.miss <- sum(is.na(S2.raw[[iExon]]))

    valid <- (I1.raw[[iExon]] + S1.raw[[iExon]] != 0) & (I2.raw[[iExon]] + S2.raw[[iExon]] !=0)
    valid[is.na(valid)] <- FALSE
      
    if(I1.miss == M[iExon] | S1.miss == M[iExon] | I2.miss == M[iExon] | S2.miss == M[iExon] | sum(valid) == 0){
      miss_total[missIndex] <- iExon
      missIndex <- missIndex + 1
    }
  }
  
  ## If all the values are missing, skip this exon/event.
  if(length(miss_total) > 0){
    I1.raw <- I1.raw[-miss_total]
    S1.raw <- S1.raw[-miss_total]

    I2.raw <- I2.raw[-miss_total]
    S2.raw <- S2.raw[-miss_total]
  
    length_I.raw <- length_I.raw[-miss_total]
    length_S.raw <- length_S.raw[-miss_total]
  
    exonList.raw <- exonList.raw[-miss_total]
    
    M <- M[-miss_total]; 

  }
  
  nExon <- length(I1.raw)
  output <- list(I1.raw, S1.raw, I2.raw, S2.raw, length_I.raw, length_S.raw, exonList.raw, nExon, M)
  names(output) <- c("I1", "S1", "I2", "S2", "length_I", "length_S", "exonList", "nExon", "M")
  
  output
  
}

