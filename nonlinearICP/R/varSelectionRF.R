varSelectionRF <- function(X, Y, env, verbose,
                           nSelect = sqrt(ncol(X)),
                           useMtry = sqrt(ncol(X)),
                           ntree = 100,
                           sampsize = min(nrow(X), 500)
                           ){
  p <- ncol(X)
  rfResultOnlyX <- randomForest(x = X, y = Y, mtry = useMtry,
                                ntree = ntree, sampsize = sampsize)
  if(nSelect < p){
    idxSelect <- order(rfResultOnlyX$importance, decreasing = TRUE)[1:nSelect]

    if(verbose) cat(paste("\n\n  ----- Choosing", nSelect, "variables out of", p,
                          "variables. Chosen variables are",
                          paste(idxSelect, collapse = ", ")))
  }else{
    idxSelect <- 1:p
  }

  idxSelect
}
