getblanketall <- function(usevar, maxSize){

  testsets <- list()

  if(length(usevar)>0){
    for (ic in ((1:2^length(usevar))-1)){
      testsets[[ic+1]] <- usevar[which( ((ic %/% 2^(0:(length(usevar)-1))) %% 2 )==1)]
    }
  }

  testsets <- unique(testsets)
  le <- sapply(testsets,length)
  testsets <- testsets[order(le)]

  # exclude empty set
  testsets[[1]] <- NULL

  # exclude sets larger than maxsize
  if(maxSize < length(usevar))
    testsets <- testsets[-which(sapply(testsets, length) > maxSize)]

  return(testsets)
}
