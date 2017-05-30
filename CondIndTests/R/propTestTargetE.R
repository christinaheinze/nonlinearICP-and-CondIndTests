#' Proportion test to compare two misclassification rates.
#'
#' @description Used as a subroutine in \code{InvariantEnvironmentPrediction} to test
#' whether out-of-sample performance is better when using X and Y as predictors for E,
#' compared to using X only.
#'
#' @param E An n-dimensional vector.
#' @param predictedOnlyX Predictions for E based on predictors in X only.
#' @param predictedXY Predictions for E based on predictors in X and Y.
#' @param verbose Set to \code{TRUE} if output should be printed.
#'
#' @return A list with the p-value for the test.
propTestTargetE <- function(E, predictedOnlyX, predictedXY, verbose){

  n <- NROW(E)

  if(!is.factor(E)){
    stop("propTestTargetE can only be applied if E is a factor.")
  }

  if(NCOL(E) > 1){
    stop("propTestTargetE can only be applied if E is univariate.")
  }

  accOnlyX <- sum(predictedOnlyX == E)
  accXY <- sum(predictedXY == E)
  testResult <- prop.test(c(accXY, accOnlyX), c(n, n), alternative="greater")
  pvalue <- testResult$p.value

  if(verbose)
    cat(paste("\nTest statistc: ", testResult$statistic))

  if(verbose)
    cat(paste("\np-value: ", pvalue))

  list(pvalue = pvalue)
}
