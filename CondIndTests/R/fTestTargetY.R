#' F-test for a nested model comparison.
#'
#' @description Used as a subroutine in \code{InvariantTargetPrediction} to test
#' whether out-of-sample prediction performance is better when using X and E as predictors for Y,
#' compared to using X only.
#'
#' @param Y An n-dimensional vector.
#' @param predictedOnlyX Predictions for Y based on predictors in X only.
#' @param predictedXE Predictions for Y based on predictors in X and E.
#' @param verbose Set to \code{TRUE} if output should be printed.
#' @param ... The dimensions of X (df) and E (dimE) need to be passed via the ...
#' argument to allow for coherent interface of fTestTargetY and wilcoxTestTargetY.
#'
#' @return A list with the p-value for the test.
fTestTargetY <- function(Y, predictedOnlyX, predictedXE, verbose, ...){
  n <- NROW(Y)
  dots <- list(...)
  if(!is.null(dots$dimE)){
    dimE <- dots$dimE
    p <- dots$df
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

  if(verbose)
    cat(paste("\np-value: ", pvalue))

  list(pvalue = pvalue)
}
