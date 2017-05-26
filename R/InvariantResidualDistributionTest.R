#' Invariant residual distribution test.
#'
#' @description Tests the null hypothesis that Y and E are independent given X.
#'
#' @param Y An n-dimensional vector.
#' @param E An n-dimensional factor variable.
#' @param X A matrix or dataframe with n rows and p columns.
#' @param alpha Significance level. Defaults to 0.05.
#' @param verbose If \code{TRUE}, intermediate output is provided. Defaults to \code{FALSE}.
#' @param fitWithGam If \code{TRUE}, a GAM is used for the nonlinear regression, else
#' a random forest is used.
#' @param test Unconditional independence test that tests whether residual distribution is
#' invariant across different levels of E. Defaults to \code{leveneAndWilcoxResidDistributions}.
#' @param colNameNoSmooth Gam parameter: Name of variables that should enter linearly into the model.
#' @param mtry Random forest parameter: Number of variables randomly sampled as
#' candidates at each split. Defaults to sqrt(NCOL(X)).
#' @param ntree Random forest parameter: Number of trees to grow. Defaults to 500.
#' @param nodesize Random forest parameter: Minimum size of terminal nodes. Defaults to 5.
#' @param maxnodes Random forest parameter: Maximum number of terminal nodes trees in the forest can have.
#' Defaults to NULL.
#' @param nSeqTests Bonferroni adjustment factor if previous tests where performed
#' (e.g. with subsamples).
#' @param returnModel If \code{TRUE}, the fitted quantile regression forest model
#' will be returned.
#'
#' @return A list with the following entries:
#' \itemize{
#'  \item \code{pValue} The p value for the null hypothesis that Y and E are independent given X.
#'  \item \code{model} The fitted model if \code{returnModel = TRUE}.
#'  }
#'
#' @examples
#' n <- 1000
#' E <- rbinom(n, size = 1, prob = 0.2)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * (X)^2 + rnorm(n)
#' InvariantResidualDistributionTest(Y, as.factor(E), X)
#' InvariantResidualDistributionTest(Y, as.factor(E), X, test = ksErrorDistributions)

#' E <- rbinom(n, size = 1, prob = 0.2)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * E + rnorm(n)
#' InvariantResidualDistributionTest(Y, as.factor(E), X)
#' InvariantResidualDistributionTest(Y, as.factor(E), X, test = ksErrorDistributions)

InvariantResidualDistributionTest <- function(Y, E, X,
                                              alpha = 0.05,
                                              verbose = FALSE,
                                              fitWithGam = TRUE,
                                              test = leveneAndWilcoxResidDistributions,
                                              colNameNoSmooth = NULL,
                                              mtry = sqrt(NCOL(X)),
                                              ntree = 500,
                                              nodesize = 5,
                                              maxnodes = NULL,
                                              nSeqTests = 1,
                                              returnModel = FALSE){

  if(!is.factor(E)){
    stop("InvariantResidualDistributionTest can only be applied if E is a factor.")
  }

  if(NCOL(E) > 1){
    stop("InvariantResidualDistributionTest can only be applied if E is univariate.")
  }

  n <- nrow(X)
  p <- ncol(X)

  if(fitWithGam){
    res <- gamErrorDistributions(X, Y, colNameNoSmooth, returnModel)
  }else{
    res <- rfErrorDistributions(X, Y, mtry, ntree, nodesize, maxnodes, returnModel)
  }

  # test whether residual distribution is identical in all environments E
  result <- test(Y, res$predicted, E, n, p, alpha, nSeqTests, verbose)


  if(returnModel){
    result$model <- res$model
  }

  result
}

