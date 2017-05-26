getgamPredictions <- function(X, Y, XTest, idxNoSmooth = NULL, returnModel = FALSE){

  p <- ncol(X)
  colnames(X) <- colnames(XTest) <- paste("V",as.character(1:p),sep="")
  cX <- colnames(X)
  data <- data.frame(Y, X)

  tr <- try({
    if(is.null(idxNoSmooth)){
      if(ncol(X) == 1 & all(X == 1)){
        form <- as.formula("Y ~ 1")
      }else{
        form <- as.formula( paste( "Y ~ ", paste( "s(", cX,")",collapse=" + "),sep=""))
      }

    }else{
      if(ncol(X) == 2 & all(X[,1] == 1)){
        form <- as.formula( paste( c("Y ~ ", paste(cX[idxNoSmooth],collapse=" + ")),collapse="") )
      }else{
        form <- as.formula( paste( c("Y ~ ", paste("s(", cX[setdiff(1:p,idxNoSmooth)],")",
                                                 collapse=" + "), " + " ,
                                   paste(cX[idxNoSmooth],collapse=" + ")),collapse="") )
      }
    }

    ga <- mgcv::gam( form, data=data, family="gaussian")

  })

  if(class(tr)[1]=="try-error"){
    form <- as.formula(Y ~ .)
    cat("\n *** reducing smoothing df")
    ga <- lm(form, data=data)
  }

  predictions <- predict(ga, newdata = XTest)
  list(predictions = predictions, model = if(returnModel) ga else NULL)
}
