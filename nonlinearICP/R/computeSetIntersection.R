#' Function to compute the intersection of sets, given a list with sets of indices.
#'
#' @param acceptedSets List containing sets of indices.
#' @param verbose Boolean variable to indicate whether messages should be printed.
#'
#' @return Vector with indices contained in the intersection of the sets.
computeSetIntersection <- function(acceptedSets, verbose = FALSE){

  if(length(acceptedSets) >= 2){

    for(i in 2:length(acceptedSets)){
      if(i == 2) current.set <- acceptedSets[[i-1]]

      current.set <- intersect(current.set, acceptedSets[[i]])
      if(length(current.set) == 0){
        if(verbose) cat("\n\n  ----- Intersection is empty, returning empty set  ----- ")
        break
      }
    }
  }else if(length(acceptedSets) == 1){
    current.set <- acceptedSets[[1]]
  }else{
    # if acceptedSets is empty:
    current.set <- numeric(0)
  }

  current.set
}
