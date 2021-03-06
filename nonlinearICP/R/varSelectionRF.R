#' Variable selection function that can be provided to \code{nonlinearICP} - it
#' is then applied to pre-select a set of variables before running the ICP procedure
#' on this subset. Here, the variable selection is based on random forest variable
#' importance measures.
#'
#' @param X A (nxp)-dimensional matrix (or data frame) with n observations of p variables.
#' @param Y Response vector (n x 1)
#' @param env Indicator of the experiment or the intervention type an observation belongs to.
#' A numeric vector of length n. Has to contain at least two different unique values.
#' @param verbose If \code{FALSE}, most messages are supressed.
#' @param nSelect Number of variables to select. Defaults to \code{sqrt(ncol(X))}.
#' @param useMtry Random forest parameter \code{mtry}. Defaults to \code{sqrt(ncol(X))}.
#' @param ntree Random forest parameter \code{ntree}. Defaults to 100.
#'
#' @return A vector containing the indices of the selected variables.
#'
#' @examples
#' # Example 1
#' require(CondIndTests)
#' data("simData")
#' targetVar <- 2
#' # choose environments where we did not intervene on var
#' useEnvs <- which(simData$interventionVar[,targetVar] == 0)
#' ind <- is.element(simData$environment, useEnvs)
#' X <- simData$X[ind,-targetVar]
#' Y <- simData$X[ind,targetVar]
#' E <- as.factor(simData$environment[ind])
#' chosenIdx <- varSelectionRF(X = X, Y = Y, env = E, verbose = TRUE)
#' cat(paste("Variable(s)", paste(chosenIdx, collapse=", "), "was/were chosen."))
varSelectionRF <- function(X, Y, env, verbose,
                           nSelect = sqrt(ncol(X)),
                           useMtry = sqrt(ncol(X)),
                           ntree = 100){
  nSelect <- round(nSelect,0)
  p <- ncol(X)
  rfResultOnlyX <- randomForest(x = X, y = Y, mtry = useMtry, ntree = ntree)

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
