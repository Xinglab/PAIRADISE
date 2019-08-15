#' pairadise
#'
#'
#' Primary function of the PAIRADISE package. Analyzes matched pairs for differences in isoform expression.
#' Uses parallel processing to speed up computation.
#'
#' @name pairadise
#' @param my.data Data frame containing grouped data to be analyzed
#' @param numCluster Number of clusters to use for parallel computing.
#' @param sig.level Positive number between 0 and 1. Specifies the desired significance level.
#' Default value is sig.level = 0.01
#' @param nIter Positive integer. Specifies the maximum number of iterations of the optimization
#' algorithm allowed. Default is nIter = 100
#' @param tol Positive number. Specifies the tolerance level for terminating the optimization
#' algorithm, defined as the difference in log-likelihood ratios between iterations. Default
#' is tol = 10^(-2)
#' @param pseudocount Positive number. Specifies a value for a pseudocount added to each
#' count at the beginning of the analysis. Default is pseudocount = 0
#' @param equal.variance Are the group variances assumed equal? Default value is FALSE.
#' @details This is the primary function of the PAIRADISE package that implements the PAIRADISE algorithm.
#' The data frame in my.data should have 7 columns, arranged as follows:
#' Column 1 contains the ID of the exons/events.
#' Column 2 contains counts of isoform 1 corresponding to the first group.
#' Column 3 contains counts of isoform 2 corresponding to the first group.
#' Column 4 contains counts of isoform 1 corresponding to the second group.
#' Column 5 contains counts of isoform 2 corresponding to the second group.
#' Replicates in columns 2-5 should be separated by commas, e.g. 1623,432,6 for three replicates.
#' Column 6 contains the effective length of isoform 1.
#' Column 7 contains the effective length of isoform 2.
#'
#' @return The function pairadise returns a list containing the following entries:
#' \item{sig.results.Bonferroni}{Matrix containing the significant exons (after Bonferroni correction at sig.level), their p-values, and test-statistics.}
#' \item{sig.results.FDR}{Matrix containing the significant exons (after FDR correction using BH at sig.level), their p-values, and test-statistics.}
#' \item{testStats}{Vector of test statistics for paired analysis.}
#' \item{raw.pvalues}{Vector of pvalues for each exon/event.}
#' \item{param.unconstrained}{List of parameter estimates for unconstrained model.}
#' \item{param.constrained}{List of parameter estimates for constrained model.}
#' \item{latent.u}{List of parameter estimates of latent variables for unconstrained model.}
#' \item{latent.c}{List of parameter estimates of latent variables for constrained model.}
#' \item{nReplicates}{Vector containing the number of valid replicates for each exon in my.data.}
#' \item{totalIter}{Vector containing total number of iterations required for optimization algorithm.}
#' \item{exonID}{Character vector containing exonIDs .}
#' \item{nExon}{Total number of valid exons in my.data.}
#' \item{I1}{List containing all exon inclusion counts for group 1.}
#' \item{S1}{List containing all exon skipping counts for group 1.}
#' \item{I2}{List containing all exon inclusion counts for group 2.}
#' \item{S2}{List containing all exon skipping counts for group 2.}
#' @examples
#'
#' #############################
#' ## Example: Simulated data ##
#' #############################
#'
#'set.seed(12345)
#'nExon <- 3  # number of exons
#'
#'## Organize data into the data frame my.data following the proper formatting:
#'exonID <- paste("Exon", as.character(seq(1:nExon)))
#'my.data <- data.frame(matrix(nrow = nExon, ncol = 7))
#'my.data[,1] <- exonID
#'my.data[,2] <- c("12,3,5", "2,9,10,6,5,4", "15,17000,20,100")
#'my.data[,3] <- c("0,1,2", "0,0,4,0,3,2", "2,12,1,1")
#'my.data[,4] <- c("2,4,5", "12,13,7,7,7,8", "1,6,7,10")
#'my.data[,5] <- c("0,1,3", "0,0,0,4,3,1", "274,NA,320,5650")
#'my.data[,6] <- c(3,3,3)
#'my.data[,7] <- c(1,1,1)
#'
#' ## Store results
#' results <- pairadise(my.data, numCluster = 4, equal.variance = FALSE)
#'
#' @export


pairadise <- function(my.data,
                           numCluster = 2,
                           sig.level = 0.01,
                           nIter = 100,
                           tol = 10^(-2),
                           pseudocount = 0,
                           seed = 12321,
                           equal.variance = FALSE){

  if(numCluster <= 0 | numCluster > 8){
    stop("Error: Number of clusters must be between 1 and 8")
  }

  if(sig.level <= 0 | sig.level >= 1){
    stop("Error: Significance level must be strictly between 0 and 1")
  }

  if(nIter <= 0){
    stop("Error: Number of iterations must be at least 1")
  }

  if(tol <= 0){
    stop("Error: Tolerance must be strictly positive")
  }

  if(pseudocount < 0){
    stop("Error: Pseudocount must be nonnegative")
  }

  if(!(equal.variance == FALSE | equal.variance == TRUE)){
    stop("Error: equal.variance must either be TRUE or FALSE")
  }

  ## Import data
  data.list <- clean.data(my.data)
  I1.main   <- data.list$I1
  S1.main   <- data.list$S1
  I2.main   <- data.list$I2
  S2.main   <- data.list$S2
  length_I  <- data.list$length_I
  length_S  <- data.list$length_S
  exonList  <- as.character(data.list$exonList)
  nExon     <- data.list$nExon

  ## Clear space in memory
  rm(data.list)

  cat("A total of", nExon, "exons will be tested.\n")
  cat("Preparing", numCluster, "clusters for parallel processing....\n")
  cl <- makeCluster(numCluster)
  registerDoParallel(cl)
  cat("Starting analysis.\n")

  writeLines(c(""), "pairadise_status.txt")

  results <- foreach(iExon = 1:nExon, .combine = cbind, .packages = c("nloptr")) %dopar% {

    ## Write progress to log
    sink("pairadise_status.txt", append=TRUE)
    cat("ExonID ", exonList[iExon], '.[',
          paste0(round(iExon / nExon * 100), '% completed]', '\n'), sep = "")
    sink()

    l.iI <- length_I[iExon]
    l.iS <- length_S[iExon]

    valid <- ((I1.main[[iExon]] + S1.main[[iExon]] != 0) & (I2.main[[iExon]] + S2.main[[iExon]] !=0))
    valid[is.na(valid)] <- FALSE

    I1 <- I1.main[[iExon]][valid] + pseudocount
    S1 <- S1.main[[iExon]][valid] + pseudocount
    I2 <- I2.main[[iExon]][valid] + pseudocount
    S2 <- S2.main[[iExon]][valid] + pseudocount

    M <- sum(valid)

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## ~~~~~~~~~~~~~~~~#~ Step 1: Initialize parameters ~#~~~~~~~~~~~~~~~~ ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

    ## Set seed
    set.seed(seed)

    ## The suffix .c always stands for "constrained"
    ## The suffix .u always stands for "unconstrained"
    logit.psi1.u <- matrix(0, nrow = (nIter + 1), ncol = M)
    logit.psi2.u <- matrix(0, nrow = (nIter + 1), ncol = M)
    alpha.u <- matrix(0, nrow = (nIter + 1), ncol = M)

    logit.psi1.c <- matrix(0, nrow = (nIter + 1), ncol = M)
    logit.psi2.c <- matrix(0, nrow = (nIter + 1), ncol = M)
    alpha.c <- matrix(0, nrow = (nIter + 1), ncol = M)

    ## eps to prevent logit(0) or logit(1)
    eps <- 0.05
    logit.psi1.u[1,] <- logit((I1 * l.iS + eps) / (I1 * l.iS + S1 * l.iI + 2 * eps))
    logit.psi2.u[1,] <- logit((I2 * l.iS + eps) / (I2 * l.iS + S2 * l.iI + 2 * eps))
    alpha.u[1,] <- logit.psi1.u[1,] + rnorm(M, mean = 0, sd = 0.01)

    logit.psi1.c[1,] <- logit.psi1.u[1,]
    logit.psi2.c[1,] <- logit.psi2.u[1,]
    alpha.c[1,] <- alpha.u[1,]

    s1.u <- c(); s1.c <- c()
    s2.u <- c(); s2.c <- c()
    s.u <- c();  s.c <- c()
    mu.u <- c(); mu.c <- c()
    delta.u <- c(); delta.c <- c()

    s1.u[1] <- 0.1;  s1.c[1] <- 0.1
    s2.u[1] <- 0.1;  s2.c[1] <- 0.1
    s.u[1]  <- 0.1;   s.c[1] <- 0.1
    mu.u[1] <- 0;    mu.c[1] <- 0
    delta.u[1] <- 0; delta.c[1] <- 0

    logit.psi1.u.old <- logit.psi1.u[1,]
    logit.psi2.u.old <- logit.psi2.u[1,]
    alpha.u.old <- alpha.u[1,]

    logit.psi1.c.old <- logit.psi1.c[1,]
    logit.psi2.c.old <- logit.psi2.c[1,]
    alpha.c.old <- alpha.c[1,]

    s1.u.old <- s1.u[1];       s1.c.old <- s1.c[1]
    s2.u.old <- s2.u[1];       s2.c.old <- s2.c[1]
    s.u.old <- s.u[1];         s.c.old <- s.c[1]
    mu.u.old <- mu.u[1];       mu.c.old <- mu.c[1]
    delta.u.old <- delta.u[1]; delta.c.old <- delta.c[1]

    ll.old.u <- 10^(40)
    ll.old.c <- 10^(40)

    ## Define limits of parameters
    s.lower <- 10^(-6); s.upper <- Inf
    mu.lower <- -Inf; mu.upper <- Inf

    MLE1.lower.u <- c(s.lower, s.lower, s.lower, mu.lower, -Inf)
    MLE1.lower.c <- c(s.lower, s.lower, s.lower, mu.lower, 0)

    MLE1.upper.u <- c(s.upper, s.upper, s.upper, mu.upper, Inf)
    MLE1.upper.c <- c(s.upper, s.upper, s.upper, mu.upper, 0)

    for(t in 1:nIter){

      ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
      ## ~~~~~~~~~~~#~ Step 2: Estimate delta, mu, and sigmas ~#~~~~~~~~~~~~ ##
      ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

      ## This is the second step in the optimization process where we estimate
      ## the MLEs  of delta, mu, sigma, sigma1, sigma2 based on the MLEs of
      ## logit(psi1), logit(psi2), and alpha computed in the previous stage.
      ## See optimize1 for more details

      ## .u's always referred to unconstrained MLE
      ## .c's always referred to constrained MLE

      MLE1.u <- function(x){

       optimize1(x, M, I1, S1, I2, S2, l.iI, l.iS, logit.psi1.u[t,], logit.psi2.u[t,], alpha.u[t,], equal.variance)

      }

      MLE1.c <- function(x){

       optimize1(x, M, I1, S1, I2, S2, l.iI, l.iS, logit.psi1.c[t,], logit.psi2.c[t,], alpha.c[t,], equal.variance)

      }

      ## Make sure parameters fall within optimization bounds
      s1.u.old <- max(s.lower, s1.u.old)
      s2.u.old <- max(s.lower, s2.u.old)
      s.u.old  <- max(s.lower, s.u.old)

      s1.c.old <- max(s.lower, s1.c.old)
      s2.c.old <- max(s.lower, s2.c.old)
      s.c.old  <- max(s.lower, s.c.old)

      if(equal.variance == FALSE){
        res1.u <- bobyqa(c(s1.u.old, s2.u.old, s.u.old, mu.u.old, delta.u.old), MLE1.u, lower = MLE1.lower.u, upper = MLE1.upper.u)
        res1.c <- bobyqa(c(s1.c.old, s2.c.old, s.c.old, mu.c.old, 0), MLE1.c, lower = MLE1.lower.c, upper = MLE1.upper.c)

        s1.u[t]    <- res1.u$par[1]
        s2.u[t]    <- res1.u$par[2]
        s.u[t]     <- res1.u$par[3]
        mu.u[t]    <- res1.u$par[4]
        delta.u[t] <- res1.u$par[5]

        s1.c[t]    <- res1.c$par[1]
        s2.c[t]    <- res1.c$par[2]
        s.c[t]     <- res1.c$par[3]
        mu.c[t]    <- res1.c$par[4]
        delta.c[t] <- res1.c$par[5]

      }else{
        res1.u <- bobyqa(c(s1.u.old, s.u.old, mu.u.old, delta.u.old), MLE1.u, lower = MLE1.lower.u[c(1,3,4,5)], upper = MLE1.upper.u[c(1,3,4,5)])
        res1.c <- bobyqa(c(s1.c.old, s.c.old, mu.c.old, 0), MLE1.c, lower = MLE1.lower.c[c(1,3,4,5)], upper = MLE1.upper.c[c(1,3,4,5)])

        s1.u[t]    <- res1.u$par[1]
        s2.u[t]    <- res1.u$par[1]
        s.u[t]     <- res1.u$par[2]
        mu.u[t]    <- res1.u$par[3]
        delta.u[t] <- res1.u$par[4]

        s1.c[t]    <- res1.c$par[1]
        s2.c[t]    <- res1.c$par[1]
        s.c[t]     <- res1.c$par[2]
        mu.c[t]    <- res1.c$par[3]
        delta.c[t] <- res1.c$par[4]
      }


      ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
      ## ~~~~~~~~~~~#~ Step 3: Estimate logit psi's and alpha ~#~~~~~~~~~~~~ ##
      ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

      #######################
      # x[1] = logit_psi1
      # x[2] = logit_psi2
      # x[3] = alpha
      #######################

      ## This is the third step in the optimization process where we estimate
      ## the MLEs of logit(psi1), logit(psi2), and alpha, based on the MLEs of
      ## delta, mu, sigma, sigma1, sigma2 computed in the previous stage.
      ## See optimize2 for more details

      ## .u's always referred to unconstrained MLE
      ## .c's always referred to constrained MLE

      for(k in 1:M){

        MLE2.u <- function(x){

          optimize2(x, k, I1, S1, I2, S2, l.iI, l.iS, delta.u[t], mu.u[t], s1.u[t], s2.u[t], s.u[t])

        }

        MLE2.c <- function(x){

          optimize2(x, k, I1, S1, I2, S2, l.iI, l.iS, delta.c[t], mu.c[t], s1.c[t], s2.c[t], s.c[t])

        }

        res2.u <- optim(c(logit.psi1.u.old[k], logit.psi2.u.old[k], alpha.u.old[k]), MLE2.u)
        res2.c <- optim(c(logit.psi1.c.old[k], logit.psi2.c.old[k], alpha.c.old[k]), MLE2.c)

        logit.psi1.u[t+1,k] <- res2.u$par[1]
        logit.psi2.u[t+1,k] <- res2.u$par[2]
        alpha.u[t+1,k]      <- res2.u$par[3]

        logit.psi1.c[t+1,k] <- res2.c$par[1]
        logit.psi2.c[t+1,k] <- res2.c$par[2]
        alpha.c[t+1,k]      <- res2.c$par[3]

      }


      ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
      ## ~~~~~~~~~~~~~~~#~ Step 4: Evaluate log likelihood ~#~~~~~~~~~~~~~~~ ##
      ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

      ll.new.u <- loglikelihood(M, I1, S1, I2, S2, l.iI, l.iS, logit.psi1.u[t,], logit.psi2.u[t,], alpha.u[t,], s1.u[t], s2.u[t], s.u[t], mu.u[t], delta.u[t])
      ll.new.c <- loglikelihood(M, I1, S1, I2, S2, l.iI, l.iS, logit.psi1.c[t,], logit.psi2.c[t,], alpha.c[t,], s1.c[t], s2.c[t], s.c[t], mu.c[t], delta.c[t])

      if((abs(ll.new.u - ll.old.u) < tol) & (abs(ll.new.c - ll.old.c) < tol)){
        break
      }

      ll.old.u <- ll.new.u
      ll.old.c <- ll.new.c

      logit.psi1.u.old <- logit.psi1.u[t+1,]
      logit.psi2.u.old <- logit.psi2.u[t+1,]
      alpha.u.old      <- alpha.u[t+1,]

      logit.psi1.c.old <- logit.psi1.c[t+1,]
      logit.psi2.c.old <- logit.psi2.c[t+1,]
      alpha.c.old      <- alpha.c[t+1,]

      s1.u.old    <- s1.u[t]
      s2.u.old    <- s2.u[t]
      s.u.old     <- s.u[t]
      mu.u.old    <- mu.u[t]
      delta.u.old <- delta.u[t]

      s1.c.old    <- s1.c[t]
      s2.c.old    <- s2.c[t]
      s.c.old     <- s.c[t]
      mu.c.old    <- mu.c[t]
      delta.c.old <- delta.c[t]

      }

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## ~~~~~~~~~~~~~~~~#~ Step 5: Save all values ~#~~~~~~~~~~~~~~~ ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

    ## Sometimes test.statistics can be negative due to numerical reasons.
    ## In these cases, set them to a very small number, so they correspond
    ## to p-value = 1.
    test.stat <- -2 * (ll.new.c - ll.new.u)
    test.stat[test.stat < 0] <- 10^(-6)

    output <- list(test.stat, delta.u[t], mu.u[t], mu.c[t], s1.u[t], s2.u[t], s.u[t],
                   s1.c[t], s2.c[t], s.c[t], sigmoid(logit.psi1.u[t+1,]),
                   sigmoid(logit.psi2.u[t+1,]), alpha.u[t+1,], sigmoid(logit.psi1.c[t+1,]),
                   sigmoid(logit.psi2.c[t+1,]), alpha.c[t+1,], M, t, exonList[iExon],
                   as.vector(I1 - pseudocount), as.vector(S1 - pseudocount),
                   as.vector(I2 - pseudocount), as.vector(S2 - pseudocount))

    output

    }

  stopCluster(cl)
  cat("Clusters closed\n")

  paired.testStats <- unlist(results[1,])
  deltas    <- results[2,]
  pred.mu.u <- results[3,]
  pred.mu.c <- results[4,]
  pred.s1.u <- results[5,]
  pred.s2.u <- results[6,]
  pred.s.u  <- results[7,]
  pred.s1.c <- results[8,]
  pred.s2.c <- results[9,]
  pred.s.c  <- results[10,]
  psi1.u    <- results[11,]
  psi2.u    <- results[12,]
  alpha.u   <- results[13,]
  psi1.c    <- results[14,]
  psi2.c    <- results[15,]
  alpha.c   <- results[16,]
  num.replicates <- unlist(results[17,])
  totalIter <- unlist(results[18,])
  exonID   <- unlist(results[19,])
  I1.total <- (results[20,])
  S1.total <- (results[21,])
  I2.total <- (results[22,])
  S2.total <- (results[23,])

  ## Give names
  names(paired.testStats) <- exonID
  names(num.replicates) <- exonID
  names(totalIter) <- exonID
  names(I1.total) <- exonID
  names(S1.total) <- exonID
  names(I2.total) <- exonID
  names(S2.total) <- exonID


  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
  ## ~~~~~~~~~~~~~~~~~#~ Find significant exons ~#~~~~~~~~~~~~~~~ ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

  ## All of the p-values of paired test
  total.pvals.paired <- 1 - pchisq(paired.testStats, 1)
  names(total.pvals.paired) <- exonID

  ## Use Bonferonni adjustment to find significant exons
  paired.index <- which(paired.testStats > qchisq(1 - sig.level/nExon, 1))
  paired.sig   <- exonList[paired.index]
  paired.pvals <- 1 - pchisq(paired.testStats[paired.index], 1)

  sig.results.Bonferroni <- cbind(paired.sig, paired.pvals, paired.testStats[paired.index])
  colnames(sig.results.Bonferroni) <- c("Exon ID", "Raw p-values", "Test Statistic")
  if(dim(sig.results.Bonferroni)[1] == 0){
    sig.results.Bonferroni <- "No significant exons"
  }

  ## Use FDR adjustment (Benjamini-Hochberg method) to find significant exons
  BH.pvals <- p.adjust(total.pvals.paired, method="BH")
  sig.results.FDR <- cbind(exonList[BH.pvals < sig.level], total.pvals.paired[BH.pvals < sig.level],
                           paired.testStats[BH.pvals < sig.level])
  colnames(sig.results.FDR) <- c("Exon ID", "Raw p-values", "Test Statistic")
  if(dim(sig.results.FDR)[1] == 0){
    sig.results.FDR <- "No significant exons"
  }

  ## Save latent variables psi1, psi2, alpha
  ## Save mu, s1, s2, s, delta
  latent.u <- list()
  latent.c <- list()
  param.unconstrained <- list()
  param.constrained   <- list()
  for(iExon in 1:nExon){
    latent.u[[iExon]] <- rbind(unlist(psi1.u[iExon]), unlist(psi2.u[iExon]), unlist(alpha.u[iExon]))
    latent.c[[iExon]] <- rbind(unlist(psi1.c[iExon]), unlist(psi2.c[iExon]), unlist(alpha.c[iExon]))

    param.unconstrained[[iExon]] <- rbind(unlist(pred.mu.u[iExon]), unlist(pred.s1.u[iExon]),
                                          unlist(pred.s2.u[iExon]), unlist(pred.s.u[iExon]),
                                          unlist(deltas[iExon]))

    param.constrained[[iExon]] <- rbind(unlist(pred.mu.c[iExon]), unlist(pred.s1.c[iExon]),
                                        unlist(pred.s2.c[iExon]), unlist(pred.s.c[iExon]))

    rownames(latent.u[[iExon]]) <- c("psi1.u", "psi2.u", "alpha.u")
    rownames(latent.c[[iExon]]) <- c("psi1.c", "psi2.c", "alpha.c")
    colnames(latent.u[[iExon]]) <- NULL
    colnames(latent.c[[iExon]]) <- NULL

    rownames(param.unconstrained[[iExon]]) <- c("mu.u", "s1.u", "s2.u", "s.u", "delta")
    rownames(param.constrained[[iExon]])   <- c("mu.c", "s1.c", "s2.c", "s.c")
    colnames(param.unconstrained[[iExon]]) <- NULL
    colnames(param.constrained[[iExon]])   <- NULL

  }

  names(latent.u) <- exonID
  names(latent.c) <- exonID
  names(param.unconstrained) <- exonID
  names(param.constrained)   <- exonID

  ## Final results
  result <- list(sig.results.Bonferroni, sig.results.FDR, paired.testStats, total.pvals.paired,
                 param.unconstrained, param.constrained, latent.u, latent.c, num.replicates,
                 totalIter, exonID, nExon, I1.total, S1.total, I2.total, S2.total)

  names(result) <- c("sig.results.Bonferroni", "sig.results.FDR", "testStats", "raw.pvalues",
                     "param.unconstrained", "param.constrained", "latent.u", "latent.c",
                     "nReplicates", "totalIter", "exonID", "nExon", "I1", "S1", "I2", "S2")

  result

}



