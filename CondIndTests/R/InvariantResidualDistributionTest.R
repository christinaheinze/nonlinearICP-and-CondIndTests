#' Invariant residual distribution test.
#'
#' @description Tests the null hypothesis that Y and E are independent given X.
#'
#' @param Y An n-dimensional vector.
#' @param E An n-dimensional vector. E needs to be a factor.
#' @param X A matrix or dataframe with n rows and p columns.
#' @param alpha Significance level. Defaults to 0.05.
#' @param verbose If \code{TRUE}, intermediate output is provided. Defaults to \code{FALSE}.
#' @param fitWithGam If \code{TRUE}, a GAM is used for the nonlinear regression, else
#' a random forest is used. Defaults to \code{TRUE}.
#' @param test Unconditional independence test that tests whether residual distribution is
#' invariant across different levels of E. Defaults to \code{leveneAndWilcoxResidDistributions}.
#' @param colNameNoSmooth Gam parameter: Name of variables that should enter linearly into the model.
#' Defaults to \code{NULL}.
#' @param mtry Random forest parameter: Number of variables randomly sampled as
#' candidates at each split. Defaults to \code{sqrt(NCOL(X))}.
#' @param ntree Random forest parameter: Number of trees to grow. Defaults to 100.
#' @param nodesize Random forest parameter: Minimum size of terminal nodes. Defaults to 5.
#' @param maxnodes Random forest parameter: Maximum number of terminal nodes trees in the forest can have.
#' Defaults to \code{NULL}.
#' @param returnModel If \code{TRUE}, the fitted quantile regression forest model
#' will be returned. Defaults to \code{FALSE}.
#'
#' @return A list with the following entries:
#' \itemize{
#'  \item \code{pvalue} The p-value for the null hypothesis that Y and E are independent given X.
#'  \item \code{model} The fitted model if \code{returnModel = TRUE}.
#'  }
#'
#' @examples
#'
#' # Example 1
#' n <- 1000
#' E <- rbinom(n, size = 1, prob = 0.2)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * (X)^2 + rnorm(n)
#' InvariantResidualDistributionTest(Y, as.factor(E), X)
#' InvariantResidualDistributionTest(Y, as.factor(E), X, test = ksResidualDistributions)
#'
#' # Example 2
#' E <- rbinom(n, size = 1, prob = 0.2)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * E + rnorm(n)
#' InvariantResidualDistributionTest(Y, as.factor(E), X)
#' InvariantResidualDistributionTest(Y, as.factor(E), X, test = ksResidualDistributions)

InvariantResidualDistributionTest <- function(Y, E, X,
                                              alpha = 0.05,
                                              verbose = FALSE,
                                              fitWithGam = TRUE,
                                              test = leveneAndWilcoxResidualDistributions,
                                              colNameNoSmooth = NULL,
                                              mtry = sqrt(NCOL(X)),
                                              ntree = 100,
                                              nodesize = 5,
                                              maxnodes = NULL,
                                              returnModel = FALSE){

  
  Y <- check_input_single(Y, return_vec = TRUE, str = "Y")
  E <- check_input_single(E, check_factor = TRUE, return_vec = TRUE, str = "E")
  X <- check_input_single(X, return_vec = FALSE)
  
  if(!is.factor(E)){
    stop("InvariantResidualDistributionTest can only be applied if E is a factor.")
  }

  if(NCOL(E) > 1){
    stop("InvariantResidualDistributionTest can only be applied if E is univariate.")
  }

  n <- NROW(X)
  p <- NCOL(X)

  if(fitWithGam){
    res <- gamResidualDistributions(X, Y, colNameNoSmooth, returnModel)
  }else{
    res <- rfResidualDistributions(X, Y, mtry, ntree, nodesize, maxnodes, returnModel)
  }

  # test whether residual distribution is identical in all environments E
  result <- test(Y, res$predicted, E, verbose)


  if(returnModel){
    result$model <- res$model
  }

  result
}

