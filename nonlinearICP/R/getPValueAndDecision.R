getPValueAndDecision <- function(X, Y, environment,
                                 condIndTest = gamTargetY,
                                 argsCondIndTest = NULL,
                                 alpha=0.05,
                                 varPreSelectionFunc = varSelectionRF,
                                 argsVarPreSelectionFunc = NULL,
                                 maxSizeSets = ncol(X),
                                 condIndTestNames = NULL,
                                 speedUp = TRUE,
                                 subsampleSize = c(0.1, 0.25, 0.5, 0.75, 1),
                                 retrieveDefiningsSets = TRUE,
                                 verbose = FALSE){

  n <- length(Y)

  if(speedUp){
    for(ss in seq_along(subsampleSize)){

      subsampleInd <-
        if(subsampleSize[ss] < 1)
          caTools::sample.split(environment, SplitRatio = subsampleSize[ss])
      else
        rep(TRUE, length = n)

      if(verbose)
        cat(paste("\nwith subsample of size ", sum(subsampleInd), sep=""))

      pvalueAndDecision <- try(do.call(condIndTest,
                                       c(list(X[subsampleInd,, drop = FALSE],
                                              Y[subsampleInd],
                                              environment[subsampleInd],
                                              alpha,
                                              verbose,
                                              nSeqTests = length(subsampleSize)),
                                         argsCondIndTest)),
                               silent = if(subsampleSize[ss] < 1) TRUE else FALSE)


      # if larger subsamples are still to be tested, continue
      if(inherits(pvalueAndDecision, "try-error") & subsampleSize[ss] < 1){
        if(verbose)
          cat(paste("\nError with subsample of size ", sum(subsampleInd),
                    ". Continue with larger sample size.",
                    sep=""))
        next
      }else if(inherits(pvalueAndDecision, "try-error") & subsampleSize[ss] == 1){
        # if full sample was used, throw error
        stop(paste("Error in conditional independence test.", geterrmessage()))
      }else if(pvalueAndDecision$decision == 0){
        # set has been rejected based on subsample size ss, stop testing
        break
      }
    }
  }else{
    pvalueAndDecision <- do.call(condIndTest,
                                 c(list(X,
                                        Y,
                                        environment,
                                        alpha,
                                        verbose,
                                        nSeqTests = 1),
                                   argsCondIndTest))
  }
  pvalueAndDecision
}
