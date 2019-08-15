#' load.data
#'
#'
#' Loads the matched pair data to be analyzed by PAIRADISE.
#'
#' @name load.data
#' @param my.data  Data frame containing grouped data to be analyzed. 
#' @return The function load.data returns a list containing the following entries:
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

load.data <- function(my.data){
  
  cat("Loading data...\n")
  
  if(class(my.data) != "data.frame"){
    stop("Error: Data is not in the proper format (not a data frame).")
  }
    
  if(ncol(my.data) < 7 ){
    stop("Error: Data frame has less than 7 columns. Terminating program.\n")
  }
    
  if(ncol(my.data) > 7 ){
    cat("Warning: Data frame has more than 7 columns. Proceeding with the first 7 columns.\n")
  }

  I1.raw <- as.list(as.character(my.data[,2]))
  I1.raw <- suppressWarnings(lapply(I1.raw, function(x) as.numeric(unlist(strsplit(x, split = ",")))))
  S1.raw <- as.list(as.character(my.data[,3]))
  S1.raw <- suppressWarnings(lapply(S1.raw, function(x) as.numeric(unlist(strsplit(x, split = ",")))))
  
  I2.raw <- as.list(as.character(my.data[,4]))
  I2.raw <- suppressWarnings(lapply(I2.raw, function(x) as.numeric(unlist(strsplit(x, split = ",")))))
  S2.raw <- as.list(as.character(my.data[,5]))
  S2.raw <- suppressWarnings(lapply(S2.raw, function(x) as.numeric(unlist(strsplit(x, split = ",")))))
  
  length_I.raw <- as.numeric(my.data[,6])
  length_S.raw <- as.numeric(my.data[,7])
  
  M1 <- unlist(lapply(I1.raw, length))
  M2 <- unlist(lapply(S1.raw, length))
  M3 <- unlist(lapply(I2.raw, length))
  M4 <- unlist(lapply(S2.raw, length))
  
  if(identical(M1, M2) == FALSE | identical(M1, M3) == FALSE | identical(M1, M4) == FALSE | 
       identical(M2, M3) == FALSE | identical(M2, M4) == FALSE | identical(M3, M4) == FALSE){
    stop("Error with data: some data are not matched pairs.")
  }

  exonList.raw <- my.data[,1]
  nExon.raw <- nrow(my.data)
  M <- M1
  output <- list(I1.raw, S1.raw, I2.raw, S2.raw, length_I.raw, length_S.raw, exonList.raw, nExon.raw, M)
  names(output) <- c("I1.raw", "S1.raw", "I2.raw", "S2.raw", "length_I.raw", "length_S.raw", "exonList.raw", "nExon.raw", "M")
  
  output
  
}
