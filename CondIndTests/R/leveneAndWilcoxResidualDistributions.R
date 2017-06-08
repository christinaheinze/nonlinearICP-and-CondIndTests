#' Levene and wilcoxon test to compare first and second moments of residual distributions
#'
#' @description Used as a subroutine in \code{InvariantResidualDistributionTest}
#' to test whether residual distribution remains invariant across different levels
#' of E.
#'
#' @param Y An n-dimensional vector.
#' @param predicted An n-dimensional vector of predictions for Y.
#' @param E An n-dimensional vector. \code{E} needs to be a factor.
#' @param verbose Set to \code{TRUE} if output should be printed.
#'
#' @return A list with the p-value for the test.
leveneAndWilcoxResidualDistributions <- function(Y, predicted, E, verbose){
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

  pvalue <- 2*min(pvalueWilcoxon, pvalueLevene)

  if(verbose)
    cat(paste("\np-value: ", pvalue))

  list(pvalue = pvalue)
}
