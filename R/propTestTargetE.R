#'  test to compare two misclassification rates.
#'
#' @description Tests whether ...
#'
#' @param E
#' @param predictedOnlyX
#' @param predictedXY
#' @param n
#' @param p
#' @param alpha
#' @param nSeqTests
#' @param verbose
#'
#' @return A list with the p-value for the test.
propTestTargetE <- function(E, predictedOnlyX, predictedXY, n, p, alpha, nSeqTests, verbose){

  if(!is.factor(E)){
    stop("propTestTargetE can only be applied if E is a factor.")
  }

  if(NCOL(E) > 1){
    stop("propTestTargetE can only be applied if E is univariate.")
  }

  accOnlyX <- sum(predictedOnlyX == E)
  accXY <- sum(predictedXY == E)
  testResult <- prop.test(c(accXY, accOnlyX), c(n, n), alternative="greater")
  pvalue <- testResult$p.value*nSeqTests

  if(verbose)
    cat(paste("\nTest statistc: ", testResult$statistic))

  if(verbose)
    cat(paste("\np-value: ", pvalue))

  list(pvalue = pvalue)
}
