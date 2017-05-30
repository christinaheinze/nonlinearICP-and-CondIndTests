#' Fishers test to test whether the exceedance of the conditional quantiles
#' is independent of the categorical variable E.
#'
#' @description Used as a subroutine in \code{InvariantConditionalQuantilePrediction}
#' to test whether the exceedance of the conditional quantiles
#' is independent of the categorical variable E.
#'
#' @param Y An n-dimensional vector.
#' @param predicted A matrix with n rows. The columns contain predictions for different
#' conditional quantiles of Y|X.
#' @param E An n-dimensional vector, defining the grouping.
#' @param adjFactor Bonferroni adjustment factor for p-value if multiple tests were performed.
#' @param verbose Set to \code{TRUE} if output should be printed.
#'
#' @return A list with the p-value for the test.

fishersTestExceedance <- function(Y, predicted, E, adjFactor, verbose){

  if(!is.factor(E)){
    stop("fishersTestExceedance can only be applied if E is a factor.")
  }

  if(NCOL(E) > 1){
    stop("fishersTestExceedance can only be applied if E is univariate.")
  }

  uniqueE <- unique(E)
  numUniqueE <- length(uniqueE)
  numQuantiles <- ncol(predicted)
  residSign <- sign( sweep(predicted, 1, Y, FUN="-") )
  # avoid zeros to have more power in test (less categories)
  residSign[residSign == 0] <- sample(c(1,-1), replace = TRUE, size = sum(residSign == 0))

  pvalue <- 1

  for(q in 1:numQuantiles){
    quantileResids <- residSign[,q]
    quantileResids <- factor(quantileResids, levels = c(1,-1))
    for(e in 1:numUniqueE){
      oneVsAllFactor <- factor((E == uniqueE[e]))
      tableE <- xtabs( ~ oneVsAllFactor + quantileResids)
      pvalue <- min(pvalue, fisher.test(tableE)$p.value)
      if(numUniqueE == 2) break
    }
  }

  bonfAdjustment <- (if(numUniqueE == 2) 1 else numUniqueE)*numQuantiles
  pvalue <- pvalue*bonfAdjustment*adjFactor

  if(verbose)
    cat(paste("\np-value: ", pvalue))

  list(pvalue = pvalue)
}
