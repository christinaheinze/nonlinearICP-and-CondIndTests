rfTargetY <- function(X,
                      Y,
                      E,
                      trainInd,
                      testInd,
                      verbose,
                      mtry,
                      ntree,
                      nodesize,
                      maxnodes,
                      permute,
                      returnModel){

  if(!permute){
    # train model with X only
    matX <- as.matrix(X)[trainInd, , drop = FALSE]
    colnames(matX) <- paste("V", 1:ncol(matX), sep = "")
    rfResultOnlyX <- randomForest(x = matX,
                                  y = Y[trainInd],
                                  mtry = mtry,
                                  ntree = ntree,
                                  nodesize = nodesize,
                                  maxnode = maxnodes)

    # predict out-of-sample with X-model
    matXPred <- as.matrix(X)[testInd, , drop = FALSE]
    colnames(matXPred) <- paste("V", 1:ncol(matXPred), sep = "")
    predictedOnlyX <- predict(rfResultOnlyX, newdata = matXPred)
  }else{
    # train model with X and E
    matX <- cbind(X = as.matrix(X)[trainInd,], E = E[trainInd][sample(1:length(trainInd))])
    colnames(matX) <- paste("V", 1:ncol(matX), sep = "")
    rfResultOnlyX <- randomForest(x = matX,
                                  y = Y[trainInd],
                                  mtry = mtry,
                                  ntree = ntree,
                                  nodesize = nodesize,
                                  maxnode = maxnodes)

    # predict out-of-sample with XE-model
    matXPred <- cbind(X = as.matrix(X)[testInd,], E = E[testInd][sample(1:length(testInd))])
    colnames(matXPred) <- paste("V", 1:ncol(matXPred), sep = "")
    predictedOnlyX <- predict(rfResultOnlyX, newdata = matXPred)
  }


  # train model with X and E
  matXE <- cbind(X = as.matrix(X)[trainInd,], E = E[trainInd])
  colnames(matXE) <- paste("V", 1:ncol(matXE), sep = "")
  rfResultXE <- randomForest(x = matXE,
                             y = Y[trainInd],
                             mtry = mtry,
                             ntree = ntree,
                             nodesize = nodesize,
                             maxnode = maxnodes)

  # predict out-of-sample with XE-model
  matXEPred <- cbind(X = as.matrix(X)[testInd,], E = E[testInd])
  colnames(matXEPred) <- paste("V", 1:ncol(matXEPred), sep = "")
  predictedXE <- predict(rfResultXE, newdata = matXEPred)

  if(!returnModel){
    list(predictedOnlyX = predictedOnlyX,
         predictedXE = predictedXE)
  }else{
    list(predictedOnlyX = predictedOnlyX,
         predictedXE = predictedXE,
         model = list(rfResultOnlyX = rfResultOnlyX, rfResultXE = rfResultXE))
  }
}
