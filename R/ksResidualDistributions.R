#' KS test to compare residual distributions
#'
#' @description Tests whether ...
#'
#' @param Y
#' @param predicted
#' @param E
#' @param n
#' @param p
#' @param alpha
#' @param nSeqTests
#' @param verbose
#'
#' @return A list with the p-value for the test.
ksResidualDistributions <- function(Y, predicted, E, n, p, alpha, nSeqTests, verbose){

  uniqueE <- unique(E)
  numUniqueE <- length(uniqueE)
  residuals <- Y - predicted
  pvalue <- 1

  # 1-vs-all
  for(e in 1:numUniqueE){
    pvalue <- min(pvalue, ks.test( residuals[which(E == uniqueE[e])], residuals[which(E != uniqueE[e])] )$p.value)
    if(numUniqueE == 2) break
  }

  bonfAdjustment <- if(numUniqueE == 2) 1 else numUniqueE

  pvalue <- pvalue*bonfAdjustment*nSeqTests

  if(verbose)
    cat(paste("\np-value: ", pvalue))

  list(pvalue = pvalue)
}
