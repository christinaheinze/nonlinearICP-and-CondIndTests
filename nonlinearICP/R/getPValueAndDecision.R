getPValueAndDecision <- function(X, Y, environment,
                                 condIndTest,
                                 argsCondIndTest,
                                 alpha,
                                 varPreSelectionFunc,
                                 argsVarPreSelectionFunc,
                                 maxSizeSets,
                                 condIndTestNames,
                                 speedUp,
                                 subsampleSize,
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

      pvalue <- try(do.call(condIndTest,
                                       c(list(Y[subsampleInd],
                                              environment[subsampleInd],
                                              X[subsampleInd,, drop = FALSE],
                                              alpha,
                                              verbose
                                              ),
                                         argsCondIndTest)),
                               silent = if(subsampleSize[ss] < 1) TRUE else FALSE)
      pvalue$pvalue <- pvalue$pvalue*length(subsampleSize)

      # if larger subsamples are still to be tested, continue
      if(inherits(pvalue, "try-error") & subsampleSize[ss] < 1){
        if(verbose)
          cat(paste("\nError with subsample of size ", sum(subsampleInd),
                    ". Continue with larger sample size.",
                    sep=""))
        next
      }else if(inherits(pvalue, "try-error") & subsampleSize[ss] == 1){
        # if full sample was used, throw error
        stop(paste("Error in conditional independence test.", geterrmessage()))
      }else if(pvalue$decision == 0){
        # set has been rejected based on subsample size ss, stop testing
        break
      }
    }
  }else{
    pvalue <- do.call(condIndTest,
                                 c(list(Y,
                                        environment,
                                        X,
                                        alpha,
                                        verbose),
                                   argsCondIndTest))
  }

  # retain null hypothesis: decision is 1
  # reject null hypothesis: decision is 0
  decision <- if(alpha >= pvalue$pvalue) 0 else 1

  list(pvalue = pvalue$pvalue, decision = decision)
}
