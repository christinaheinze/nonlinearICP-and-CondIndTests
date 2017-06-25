library(nonlinearICP)
context("All supported methods")

data("simData")

targetVar <- 2
# choose environments where we did not intervene on var
useEnvs <- which(simData$interventionVar[,targetVar] == 0)
ind <- is.element(simData$environment, useEnvs)
X <- simData$X[ind,-targetVar]
Y <- simData$X[ind,targetVar]
E <- as.factor(simData$environment[ind])

# methods and tests
methods <- NULL
tests <- NULL

methods <- c(methods, list("InvariantTargetPrediction" = InvariantTargetPrediction,
                           "InvariantTargetPrediction" = InvariantTargetPrediction))

tests <- c(tests, list("wilcoxTestTargetY" = wilcoxTestTargetY,
                       "fTestTargetY" = fTestTargetY))

methods <- c(methods, list("InvariantResidualDistributionTest" = InvariantResidualDistributionTest,
                           "InvariantResidualDistributionTest" = InvariantResidualDistributionTest))

tests <- c(tests, list("ksErrorDistributions" = ksResidualDistributions,
                       "leveneAndWilcoxErrorDistributions" = leveneAndWilcoxResidualDistributions))

methods <- c(methods, list("InvariantConditionalQuantilePrediction" = InvariantConditionalQuantilePrediction))
tests <- c(tests, list("fishersTestExceedance" = fishersTestExceedance))

methods <- c(methods, list("ResidualPredictionTest" = ResidualPredictionTest))
tests <- c(tests, list("NULL" = NULL))

methods <- c(methods, list("InvariantEnvironmentPrediction" = InvariantEnvironmentPrediction))
tests <- c(tests, list("propTestTargetE" = propTestTargetE))

methods <- c(methods, list("KCI" = KCI))
tests <- c(tests, list("NULL" = NULL))


for(methodIdx in 1:length(methods)){
  test_that(paste("Checks output type for", names(methods)[methodIdx],
                  "and test ", names(tests)[methodIdx]), {

    argsCondIndTestList <- if(!is.null(tests[[methodIdx]])) list(test = tests[[methodIdx]])

    outNonLin <- nonlinearICP(X = X,
                          Y = Y,
                          environment = E,
                          condIndTest = methods[[methodIdx]],
                          argsCondIndTest = argsCondIndTestList,
                          condIndTestNames = c(names(methods)[methodIdx],
                                               names(tests)[methodIdx]),
                          verbose = FALSE,
                          retrieveDefiningsSets = FALSE,
                          stopIfEmpty = TRUE,
                          seed = 1)

    expect_is(
      outNonLin$acceptedSets
      , "list")

  }
  )
}
