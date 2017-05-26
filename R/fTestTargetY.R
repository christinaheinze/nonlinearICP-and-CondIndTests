#' F test to compare two fits...
#'
#' @description Tests whether ...
#'
#' @param Y
#' @param predictedOnlyX
#' @param predictedXE
#' @param n
#' @param p
#' @param dimE
#' @param alpha
#' @param nSeqTests
#' @param verbose
#'
#' @return A list with the p-value for the test.
fTestTargetY <- function(Y, predictedOnlyX, predictedXE, n, p, alpha, nSeqTests, verbose, ...){

  dots <- list(...)
  if(!is.null(dots$dimE)){
    dimE <- dots$dimE
  }else{
    dimE <- 1
    warning("Assuming E is univariate.")
  }

  sseOnlyX <- sum((Y - predictedOnlyX)^2)
  sseXE <- sum((Y - predictedXE)^2)

  fStatistic <- ((sseOnlyX-sseXE)/dimE)/(sseXE/(n-p))
  if(verbose)
    cat(paste("\nF-Statistc: ", fStatistic))

  pvalue <- pf(fStatistic, dimE, n-p, lower.tail = FALSE)
  pvalue <- nSeqTests*pvalue

  if(verbose)
    cat(paste("\np-value: ", pvalue))

  list(pvalue = pvalue)
}
