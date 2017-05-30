#' Wilcoxon test to compare two mean squared error rates.
#'
#' @description Used as a subroutine in \code{InvariantTargetPrediction} to test
#' whether out-of-sample performance is better when using X and E as predictors for Y,
#' compared to using X only.
#'
#' @param Y An n-dimensional vector.
#' @param predictedOnlyX Predictions for Y based on predictors in X only.
#' @param predictedXE Predictions for Y based on predictors in X and E.
#' @param verbose Set to \code{TRUE} if output should be printed.
#' @param ... Argument to allow for coherent interface of fTestTargetY and wilcoxTestTargetY.
#'
#' @return A list with the p-value for the test.
wilcoxTestTargetY <- function(Y, predictedOnlyX, predictedXE, verbose, ...){

  residOnlyX <- abs(Y - predictedOnlyX)
  residXE <- abs(Y - predictedXE)
  diffResid <- residOnlyX-residXE

  testResult <- wilcox.test(diffResid, alternative="greater")
  pvalue <- testResult$p.value

  if(verbose)
    cat(paste("\nTest statistc: ", testResult$statistic))

  if(verbose)
    cat(paste("\np-value: ", pvalue))

  list(pvalue = pvalue)
}
