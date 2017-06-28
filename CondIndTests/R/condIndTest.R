#' Wrapper function for conditional independence tests.
#'
#' @description Tests the null hypothesis that Y and E are independent given X.
#'
#' @param Y An n-dimensional vector or a matrix or dataframe with n rows and p columns.
#' @param E An n-dimensional vector or a matrix or dataframe with n rows and p columns.
#' @param X An n-dimensional vector or a matrix or dataframe with n rows and p columns.
#' @param method The conditional indepdence test to use, can be one of
#' \code{"KCI"}, \code{"InvariantConditionalQuantilePrediction"}, \code{"InvariantEnvironmentPrediction"},
#'  \code{"InvariantResidualDistributionTest"}, \code{"InvariantTargetPrediction"}, \code{"ResidualPredictionTest"}.
#' @param alpha Significance level. Defaults to 0.05.
#' @param parsMethod Named list to pass options to \code{method}.
#' @param verbose If \code{TRUE}, intermediate output is provided. Defaults to \code{FALSE}.
#'
#' @return A list with the p-value of the test (\code{pvalue}) and possibly additional
#' entries, depending on the output of the chosen conditional independence test in \code{method}.
#'
#' @references Please cite
#' C. Heinze-Deml, J. Peters and N. Meinshausen: "Invariant Causal Prediction for Nonlinear Models",
#' \href{https://arxiv.org/abs/1706.08576}{arXiv:1706.08576}
#' and the corresponding reference for the conditional independence test.
#'
#' @examples
#'
#' # Example 1
#' set.seed(1)
#' n <- 100
#' Z <- rnorm(n)
#' X <- 4 + 2 * Z + rnorm(n)
#' Y <- 3 * X^2 + Z + rnorm(n)
#' test1 <- CondIndTest(X,Y,Z, method = "KCI")
#' cat("These data come from a distribution, for which X and Y are NOT
#' cond. ind. given Z.")
#' cat(paste("The p-value of the test is: ", test1$pvalue))
#'
#' # Example 2
#' set.seed(1)
#' Z <- rnorm(n)
#' X <- 4 + 2 * Z + rnorm(n)
#' Y <- 3 + Z + rnorm(n)
#' test2 <- CondIndTest(X,Y,Z, method = "KCI")
#' cat("The data come from a distribution, for which X and Y are cond.
#' ind. given Z.")
#' cat(paste("The p-value of the test is: ", test2$pvalue))
#'
CondIndTest <- function(Y, E, X,
                          method = "KCI",
                          alpha = 0.05,
                          parsMethod = list(),
                          verbose = FALSE){

  argsSet <- list(Y = Y,
                  E = E,
                  X = X,
                  alpha = alpha,
                  verbose = verbose)

  switch(method,
         "KCI" = {
           result <- do.call(KCI, c(argsSet, parsMethod))
         },
         "InvariantConditionalQuantilePrediction" = {
           result <- do.call(InvariantConditionalQuantilePrediction, c(argsSet, parsMethod))
         },
         "InvariantEnvironmentPrediction" = {
           result <- do.call(InvariantEnvironmentPrediction, c(argsSet, parsMethod))
         },
         "InvariantResidualDistributionTest" = {
           result <- do.call(InvariantResidualDistributionTest, c(argsSet, parsMethod))
         },
         "InvariantTargetPrediction" = {
           result <- do.call(InvariantTargetPrediction, c(argsSet, parsMethod))
         },
         "ResidualPredictionTest" = {
           result <- do.call(ResidualPredictionTest, c(argsSet, parsMethod))
         },
         {
           stop(paste("Method ", method," not implemented"))
         }
  )

  return(result)
}
