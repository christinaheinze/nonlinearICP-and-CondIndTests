#' Invariant environment prediction.
#'
#' @description Tests the null hypothesis that Y and E are independent given X.
#'
#' @param Y An n-dimensional vector.
#' @param E An n-dimensional vector. If \code{test = propTestTargetE}, E needs
#' to be a factor.
#' @param X A matrix or dataframe with n rows and p columns.
#' @param alpha Significance level. Defaults to 0.05.
#' @param verbose If \code{TRUE}, intermediate output is provided. Defaults to \code{FALSE}.
#' @param trainTestSplitFunc Function to split sample. Defaults to stratified sampling
#' using \code{caTools::sample.split}, assuming E is a factor.
#' @param argsTrainTestSplitFunc Arguments for sampling splitting function.
#' @param test Unconditional independence test that tests whether the out-of-sample
#' prediction accuracy is the same when using X only vs. X and Y as predictors for E.
#' Defaults to \code{propTestTargetE}.
#' @param mtry Random forest parameter: Number of variables randomly sampled as
#' candidates at each split.  Defaults to \code{sqrt(NCOL(X))}.
#' @param ntree Random forest parameter: Number of trees to grow. Defaults to 100.
#' @param nodesize Random forest parameter: Minimum size of terminal nodes. Defaults to 5.
#' @param maxnodes Random forest parameter: Maximum number of terminal nodes trees in the forest can have.
#' Defaults to \code{NULL}.
#' @param permute Random forest parameter: If \code{TRUE}, model that would use X only
#' for predicting Y also includes a random permutation of E. Defaults to \code{TRUE}.
#' @param returnModel If \code{TRUE}, the fitted quantile regression forest model
#' will be returned. Defaults to \code{FALSE}.
#'
#' @return A list with the following entries:
#' \itemize{
#'  \item \code{pvalue} The p-value for the null hypothesis that Y and E are independent given X.
#'  \item \code{model} The fitted models if \code{returnModel = TRUE}.
#'  }
#'
#' @examples
#' # Example 1
#' n <- 1000
#' E <- rbinom(n, size = 1, prob = 0.2)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * (X)^2 + rnorm(n)
#' InvariantEnvironmentPrediction(Y, as.factor(E), X)
#'
#' # Example 2
#' E <- rbinom(n, size = 1, prob = 0.2)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * E + rnorm(n)
#' InvariantEnvironmentPrediction(Y, as.factor(E), X)
#'
#' # Example 3
#' E <- rnorm(n)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * (X)^2 + rnorm(n)
#' InvariantEnvironmentPrediction(Y, E, X, test = wilcoxTestTargetY)
#' InvariantEnvironmentPrediction(Y, X, E, test = wilcoxTestTargetY)

InvariantEnvironmentPrediction <- function(Y, E, X,
                                           alpha = 0.05,
                                           verbose = FALSE,
                                           trainTestSplitFunc = caTools::sample.split,
                                           argsTrainTestSplitFunc = list(Y = E, SplitRatio = 0.8),
                                           test = propTestTargetE,
                                           mtry = sqrt(NCOL(X)),
                                           ntree = 100,
                                           nodesize = 5,
                                           maxnodes = NULL,
                                           permute = TRUE,
                                           returnModel = FALSE){

  if(!is.factor(E) & length(unique(E)) < 5){
    warning("E has less than 5 unique values; are you sure that E is not a factor?")
  }

  n <- NROW(X)
  p <- NCOL(X)

  trainInd <- do.call(trainTestSplitFunc, argsTrainTestSplitFunc)
  testInd <- which(!trainInd)
  trainInd <- which(trainInd)

  if(!permute){
    if(NCOL(X) == 1 & all(X == 1) & is.factor(E)){
      # predict out-of-sample with intercept model: most frequent class in training
      # set
      envFreq <- as.data.frame(table(E))
      predictedEnv <- envFreq[which.max(envFreq$Freq),1]
      predictedOnlyX <- rep(predictedEnv, times = nrow(as.matrix(X)[testInd,]))

    }else{
      if(verbose){
        cat(paste("\nNon-permuted model with", NCOL(X), "vars in X"))
      }
      matX <- as.matrix(X)[trainInd,]
      colnames(matX) <- paste("V", 1:ncol(matX), sep = "")
      # train model with X only
      rfResultOnlyX <- randomForest(x = matX, y = E[trainInd],
                                    mtry = mtry, ntree = ntree,
                                    nodesize = nodesize, maxnodes = maxnodes)

      # predict out-of-sample with X-model
      matXPred <- as.matrix(X)[testInd,]
      colnames(matXPred) <- paste("V", 1:ncol(matXPred), sep = "")
      predictedOnlyX <- predict(rfResultOnlyX, newdata = matXPred)
    }
  }else{
    if(verbose){
      cat(paste("\nPermuted model with", NCOL(X), "vars in X"))
    }

    # train model with X and permuted Y
    matX <- cbind(X = as.matrix(X)[trainInd,], Y = Y[trainInd][sample(1:length(trainInd))])
    colnames(matX) <- paste("V", 1:ncol(matX), sep = "")
    rfResultOnlyX <- randomForest(x = matX,
                                  y = E[trainInd],
                                  mtry = mtry,
                                  ntree = ntree,
                                  nodesize = nodesize,
                                  maxnodes = maxnodes)
    # predict out-of-sample with XY-model
    matXPred <- cbind(X = as.matrix(X)[testInd,], Y = Y[testInd][sample(1:length(testInd))])
    colnames(matXPred) <- paste("V", 1:ncol(matXPred), sep = "")
    predictedOnlyX <- predict(rfResultOnlyX, newdata = matXPred)
  }

  if(verbose){
    cat(paste("\nTraining X and Y together..."))
  }

  # train model with X and Y
  matXY <- cbind(X = as.matrix(X)[trainInd,], Y = Y[trainInd])
  colnames(matXY) <- paste("V", 1:ncol(matXY), sep = "")
  rfResultXY <- randomForest(x = matXY, y = E[trainInd], mtry = mtry,
                             ntree = ntree,
                             nodesize = nodesize, maxnodes = maxnodes)

  # predict out-of-sample with XY-model
  matXYPred <- cbind(X = as.matrix(X)[testInd,], Y = Y[testInd])
  colnames(matXYPred) <- paste("V", 1:ncol(matXYPred), sep = "")
  predictedXY <- predict(rfResultXY, newdata = matXYPred)

  # test whether performance is statistically indistinguishable
  result <- test(E[testInd], predictedOnlyX, predictedXY, verbose)


  # reject if using X and E has significantly better accuracy than using X only
  if(verbose){
    if(is.factor(E)){
      cat(paste("\nAccuracy only X :", round( mean(predictedOnlyX == E[testInd]) , 2),
                "\nAccuracy with X and Y:", round( mean(predictedXY == E[testInd]) , 2)))
    }else{
      cat(paste("\nMSE only X :", round(mean((E[testInd] - predictedOnlyX)^2), 2),
                "\nMSE with X and Y:", round(mean((E[testInd] - predictedXY)^2), 2)))
    }
  }


  if(returnModel){
    result$model <- list(rfResultOnlyX = if(exists("rfResultOnlyX")) rfResultOnlyX else NULL,
                         rfResultXY = rfResultXY,
                         predictedOnlyX = predictedOnlyX,
                         predictedXY = predictedXY)
  }

  result
}
