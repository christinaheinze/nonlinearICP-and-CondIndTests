rfResidualDistributions <- function(X, Y,
                                 mtry,
                                 ntree,
                                 nodesize,
                                 maxnodes,
                                 returnModel){

  n <- nrow(X)
  p <- ncol(X)

  # train model Y ~ X using all data
  rfResult <- randomForest(x = as.matrix(X),
                           y = Y,
                           mtry = mtry,
                           ntree = ntree,
                           nodesize = nodesize,
                           maxnodes = maxnodes)

  # predict
  predicted <- predict(rfResult, newdata = as.matrix(X))

  if(!returnModel){
    list(predicted = predicted)
  }else{
    list(predicted = predicted, model = list(rfResult = rfResult))
  }
}

