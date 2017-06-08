#' Invariant conditional quantile prediction.
#'
#' @description Tests the null hypothesis that Y and E are independent given X.
#'
#' @param Y An n-dimensional vector.
#' @param E An n-dimensional vector. If \code{test = fishersTestExceedance}, E needs
#' to be a factor.
#' @param X A matrix or dataframe with n rows and p columns.
#' @param alpha Significance level. Defaults to 0.05.
#' @param verbose If \code{TRUE}, intermediate output is provided. Defaults to \code{FALSE}.
#' @param test Unconditional independence test that tests whether exceedence is
#' independent of E. Defaults to \code{fishersTestExceedance}.
#' @param mtry Random forest parameter: Number of variables randomly sampled as
#' candidates at each split. Defaults to \code{sqrt(NCOL(X))}.
#' @param ntree Random forest parameter: Number of trees to grow. Defaults to 100.
#' @param nodesize Random forest parameter: Minimum size of terminal nodes.  Defaults to 5.
#' @param maxnodes Random forest parameter: Maximum number of terminal nodes trees in the forest can have.
#' Defaults to NULL.
#' @param quantiles Quantiles for which to test independence between exceedence and E.
#' Defaults to \code{c(0.1, 0.5, 0.9)}.
#' @param returnModel If \code{TRUE}, the fitted quantile regression forest model
#' will be returned. Defaults to \code{FALSE}.
#'
#' @return A list with the following entries:
#' \itemize{
#'  \item \code{pvalue} The p-value for the null hypothesis that Y and E are independent given X.
#'  \item \code{model} The fitted quantile regression forest model if \code{returnModel = TRUE}.
#'  }
#'
#' @examples
#' # Example 1
#' n <- 1000
#' E <- rbinom(n, size = 1, prob = 0.2)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * (X)^2 + rnorm(n)
#' InvariantConditionalQuantilePrediction(Y, as.factor(E), X)
#'
#' # Example 2
#' E <- rbinom(n, size = 1, prob = 0.2)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * E + rnorm(n)
#' InvariantConditionalQuantilePrediction(Y, as.factor(E), X)
#'
InvariantConditionalQuantilePrediction <- function(Y, E, X,
                                                  alpha = 0.05,
                                                  verbose = FALSE,
                                                  test = fishersTestExceedance,
                                                  mtry = sqrt(NCOL(X)),
                                                  ntree = 100,
                                                  nodesize = 5,
                                                  maxnodes = NULL,
                                                  quantiles = c(0.1, 0.5, 0.9),
                                                  returnModel = FALSE){

  n <- NROW(X)
  p <- NCOL(X)

  # train model Y ~ X using all data
  mat <- as.matrix(X)
  colnames(mat) <- paste("V", 1:ncol(mat), sep = "")
  rfResult <- quantregForest(x = mat,
                             y = Y,
                             mtry = mtry,
                             ntree = ntree,
                             nodesize = nodesize,
                             maxnodes = maxnodes)

  # predict
  predicted <- predict(rfResult, newdata = mat, what = quantiles)

  # test whether residual distribution is identical in all environments E
  result <- test(Y, predicted, E, verbose)

  if(returnModel){
    result$model <- list(rfResult = rfResult)
  }

  result
}

