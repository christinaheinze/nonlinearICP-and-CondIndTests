wilcoxTestTargetY <- function(Y, predictedOnlyX, predictedXE, n, p, alpha, nSeqTests, verbose){

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
