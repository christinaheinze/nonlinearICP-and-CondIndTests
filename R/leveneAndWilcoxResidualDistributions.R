#'  test to compare residual distributions
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
leveneAndWilcoxResidualDistributions <- function(Y, predicted, E, n, p, alpha, nSeqTests, verbose){
  uniqueE <- unique(E)
  numUniqueE <- length(uniqueE)
  residuals <- Y - predicted
  pvalueWilcoxon <- 1

  # Levene's test for homogeneity of variance across groups
  pvalueLevene <- levene.test(residuals, as.factor(E), location = "mean")$p.value

  # Wilcoxon rank sum test for equality of expectation: 1-vs-all
  for(e in 1:numUniqueE){
    pvalueWilcoxon <- min(pvalueWilcoxon, wilcox.test( residuals[which(E == uniqueE[e])], residuals[which(E != uniqueE[e])] )$p.value)
    if(numUniqueE == 2) break
  }

  bonfAdjustment <- if(numUniqueE == 2) 1 else numUniqueE

  pvalueWilcoxon <- pvalueWilcoxon*bonfAdjustment

  pvalue <- nSeqTests*2*min(pvalueWilcoxon, pvalueLevene)

  if(verbose)
    cat(paste("\np-value: ", pvalue))

  list(pvalue = pvalue)
}
