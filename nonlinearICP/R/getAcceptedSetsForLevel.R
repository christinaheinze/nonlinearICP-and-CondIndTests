getAcceptedSetsForLevel <- function(alpha, accSets, rejSets, pvalAcc, pvalRej){

  allSets <- c(accSets, rejSets)
  allPvals <- c(pvalAcc, pvalRej)
  allSets[which(allPvals > alpha)]
}
