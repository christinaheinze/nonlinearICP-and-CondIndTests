gamResidualDistributions <- function(X, Y,
                                  colNameNoSmooth = NULL,
                                  returnModel = FALSE){

  n <- nrow(X)
  p <- ncol(X)


  if(!is.null(colNameNoSmooth)){
    idxNoSmooth <- which(is.element(colnames(X), colNameNoSmooth))

    if(length(idxNoSmooth) == 0){
      idxNoSmooth <- NULL
    }

  }else{
    idxNoSmooth <- NULL
  }

  # train model Y ~ X using all data and predict
  gamRes <- getgamPredictions(as.data.frame(as.matrix(X)),
                              Y,
                              as.data.frame(as.matrix(X)),
                              idxNoSmooth = idxNoSmooth,
                              returnModel = returnModel)
  predicted <- gamRes$predictions

  if(!returnModel){
    list(predicted = predicted)
  }else{
    list(predicted = predicted, model = gamRes$model)
  }
}

