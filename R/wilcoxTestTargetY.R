#' Wilcoxon test to compare two mean squared error rates.
#'
#' @description Tests whether ...
#'
#' @param Y
#' @param predictedOnlyX
#' @param predictedXE
#' @param n
#' @param p
#' @param alpha
#' @param nSeqTests
#' @param verbose
#'
#' @return A list with the p-value for the test.
wilcoxTestTargetY <- function(Y, predictedOnlyX, predictedXE, n, p, alpha, nSeqTests, verbose, ...){

  residOnlyX <- abs(Y - predictedOnlyX)
  residXE <- abs(Y - predictedXE)
  diffResid <- residOnlyX-residXE

  testResult <- wilcox.test(diffResid, alternative="greater")
  pvalue <- testResult$p.value*nSeqTests

  if(verbose)
    cat(paste("\nTest statistc: ", testResult$statistic))

  if(verbose)
    cat(paste("\np-value: ", pvalue))

  list(pvalue = pvalue)
}
