#' Given a significance level \code{alpha} and the results from running \code{nonlinearICP},
#' returns the accepted sets for given significance level.
#'
#' @param alpha Significance level for which to retrieve the accepted sets.
#' @param List of accepted sets as returned by \code{nonlinearICP}.
#' @param List of rejected sets as returned by \code{nonlinearICP}.
#' @param List of p-values corresponding to accepted sets as returned by \code{nonlinearICP}.
#' @param List of p-values corresponding to rejected sets as returned by \code{nonlinearICP}.
#'
#' @return List of accepted sets for significance level \code{alpha}.
getAcceptedSetsForLevel <- function(alpha, accSets, rejSets, pvalAcc, pvalRej){

  allSets <- c(accSets, rejSets)
  allPvals <- c(pvalAcc, pvalRej)
  allSets[which(allPvals > alpha)]
}
