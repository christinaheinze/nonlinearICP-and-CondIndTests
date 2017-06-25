##' Summary functions for 'nonlinICP.class' objects.
##'
##' @title summary function
##' @param object object of class 'nonlinICP.class'.
##'
##' @author Christina Heinze-Deml and Jonas Peters
##'
##' @export

summary.nonlinICP.class <- function(object, show.past=TRUE, ...){
  stopifnot(inherits(object, "nonlinICP.class"))
  numSets <- length(object$acceptedSets)
  if(numSets > 0){
    cat(paste("\n", numSets, " of the sets that were tested, were accepted."))
    cat("\nThe intersection equals:\n")
    show(object$retrievedCausalVars)
    cat("\nThis is the output of the plain method. The defining sets contain more information and might be interesting to look at --- especially if the output is empty.")
  } else {
    cat(paste("\nNone of the tests were accepted. This could be due to model misspecification. Or, the environment has a direct influence on the target variable."))
  }
}
