#' Invariant target prediction.
#'
#' @description Tests the null hypothesis that Y and E are independent given X.
#'
#' @param Y An n-dimensional vector.
#' @param E An n-dimensional vector or an nxq dimensional matrix or dataframe.
#' @param X A matrix or dataframe with n rows and p columns.
#' @param alpha Significance level. Defaults to 0.05.
#' @param verbose If \code{TRUE}, intermediate output is provided. Defaults to \code{FALSE}.
#' @param fitWithGam If \code{TRUE}, a GAM is used for the nonlinear regression, else
#' a random forest is used.
#' @param trainTestSplitFunc Function to split sample. Defaults to stratified sampling when
#' E is a factor,
#' @param argsTrainTestSplitFunc Arguments for sampling splitting function.
#' @param test Unconditional independence test that tests whether the out-of-sample
#' prediction accuracy is the same when using X only vs. X and E as predictors for Y.
#' @param colNameNoSmooth Gam parameter: Name of variables that should enter linearly into the model.
#' @param mtry Random forest parameter: Number of variables randomly sampled as
#' candidates at each split.
#' @param ntree Random forest parameter: Number of trees to grow.
#' @param nodesize Random forest parameter: Minimum size of terminal nodes.
#' @param maxnodes Random forest parameter: Maximum number of terminal nodes trees in the forest can have.
#' Defaults to NULL.
#' @param permute Random forest parameter: If \code{TRUE}, model that would use X only
#' for predicting Y also includes a random permutation of E.
#' @param nSeqTests Bonferroni adjustment factor if previous tests where performed
#' (e.g. with subsamples).
#' @param returnModel If \code{TRUE}, the fitted quantile regression forest model
#' will be returned.
#'
#' @return A list with the following entries:
#' \itemize{
#'  \item \code{pValue} The p value for the null hypothesis that Y and E are independent given X.
#'  \item \code{model} The fitted models if \code{returnModel = TRUE}.
#'  }
#'
#' @examples
#' # Example 1
#' n <- 1000
#' E <- rbinom(n, size = 1, prob = 0.2)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * (X)^2 + rnorm(n)
#' InvariantTargetPrediction(Y, as.factor(E), X)
#' InvariantTargetPrediction(Y, as.factor(E), X, test = wilcoxTestTargetY)
#'
#' # Example 2
#' E <- rbinom(n, size = 1, prob = 0.2)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * E + rnorm(n)
#' InvariantTargetPrediction(Y, as.factor(E), X)
#' InvariantTargetPrediction(Y, as.factor(E), X, test = wilcoxTestTargetY)
#'
#' # Example 3
#' E <- rnorm(n)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * (X)^2 + rnorm(n)
#' InvariantTargetPrediction(Y, E, X)
#' InvariantTargetPrediction(Y, X, E)
#' InvariantTargetPrediction(Y, E, X, test = wilcoxTestTargetY)
#' InvariantTargetPrediction(Y, X, E, test = wilcoxTestTargetY)

InvariantTargetPrediction <-  function(Y, E, X,
                                       alpha = 0.05,
                                       verbose = FALSE,
                                       fitWithGam = TRUE,
                                       trainTestSplitFunc = caTools::sample.split,
                                       argsTrainTestSplitFunc = list(SplitRatio = 0.8),
                                       test = fTestTargetY,
                                       colNameNoSmooth = NULL,
                                       mtry = sqrt(NCOL(X)),
                                       ntree = 250,
                                       nodesize = 5,
                                       maxnodes = NULL,
                                       permute = TRUE,
                                       nSeqTests = 1,
                                       returnModel = FALSE
                                       ){

  if(!is.factor(E) & length(unique(E)) < 5){
    warning("E has less than 5 unique values; are you sure that E is not a factor?")
  }

  n <- NROW(X)
  p <- NCOL(X)

  trainInd <- do.call(trainTestSplitFunc, c(list(E), argsTrainTestSplitFunc))
  testInd <- which(!trainInd)
  trainInd <- which(trainInd)

  if(verbose)
    cat(paste("\nUsing ", length(trainInd), "samples for training;", length(testInd),
              "samples for testing."))


  if(fitWithGam){
    res <- gamTargetY(X, Y, E, trainInd, testInd, verbose, colNameNoSmooth, returnModel)
  }else{
    res <- rfTargetY(X, Y, E, trainInd, testInd, verbose, mtry, ntree, nodesize, maxnodes, permute, returnModel)
  }

  # test whether performance is statistically indistinguishable
  # dfs: larger model: p + env. var. + intercept

  dimE <- NCOL(E)

  if(NCOL(X) == 1 & all(X == 1)){
    df <- 1+dimE
  }else{
    df <- p+1+dimE
  }

  result <- test(Y[testInd], res$predictedOnlyX, res$predictedXE, nSeqTests, verbose, df = df, dimE = dimE)

  # reject if using X and E has significantly better accuracy than using X only
  if(verbose) cat(paste("\nMSE only X :", round(mean((Y[testInd] - res$predictedOnlyX)^2), 2),
                        "\nMSE with X and E:", round(mean((Y[testInd] - res$predictedXE)^2), 2)))

  if(returnModel){
    result$model <- res$model
  }

  result
}

