computeDefiningSets <- function(as){

  if(length(as) == 0) return(list())

  treeList <- list()

  # add all singletons
  indSingletons <- which(sapply(as, length) == 1)
  singletons <- as[indSingletons]

  if(length(singletons) > 0){
    for(i in 1:length(singletons)){
      if(i == 1) treeList[[1]] <- Node$new(singletons[1])
      else{
        addTo <- Traverse(treeList[[1]], filterFun = isLeaf)
        addTo[[1]]$AddChild(singletons[i])
      }

      ind <- which(sapply(as, function(x) is.element(unlist(singletons[i]), x)))
      if(length(ind) > 0) as <- as[-ind]
    }

  }

  # create root node of tree(s)
  if(length(as) > 0){
    occur <- table(unlist(as))
    uniqueVars <- sort(unique(unlist(as)))

    if(length(treeList) == 0){
      for(i in 1:length(uniqueVars)){
        treeList[[length(treeList) + 1]] <- Node$new(uniqueVars[i])
      }
    }else{

      for(tl in 1:length(treeList)){
        leafNode <- Traverse(treeList[[tl]], filterFun = isLeaf)

        for(i in 1:length(uniqueVars)){
          leafNode[[1]]$AddChild(uniqueVars[i])
        }
      }
    }

    # grow trees
    for(i in 1:length(treeList)){

      cont <- TRUE

      while(cont){

        leafNodes <- Traverse(treeList[[i]], filterFun = isLeaf)
        leafNodePath <- lapply(leafNodes, function(i) i$path)
        continueAny <- rep(TRUE, times = length(leafNodePath))

        for(ln in 1:length(leafNodes)){
          # remove sets with variables already in branch
          toRemove <- as.numeric(leafNodePath[[ln]]) # leave node + path to root
          addTo <- leafNodes[[ln]]

          ind <- NULL
          for(elem in 1:length(toRemove)){
            ind <- c(ind, which(sapply(as, function(x) is.element(toRemove[elem], x))))
          }

          asReduced <- if(length(ind) > 0) as[-ind] else as

          if(length(asReduced) == 0){
            continueAny[ln] <- FALSE
          }else{
            # get subtree
            treeList[[i]] <- buildTrees(asReduced, treeList[[i]], addTo)
          }

        }
        cont <- any(continueAny)
      }

    }
  }

  # prune
  sets <- NULL
  for(i in 1:length(treeList)){
    x <- ToDataFrameTree(treeList[[i]], "pathString", "isLeaf")

    idxLeaves <- which(x$isLeaf)
    sets <- c(sets, lapply(x[idxLeaves, "pathString"],
                   function(x) sort(as.numeric(strsplit(x, "/")[[1]]))))
  }

  idxNotDup <- which(!duplicated(sets))
  setsUnique <- sets[idxNotDup]

  isSubset <- function(set1, set2){
    all(sapply(set1, function(el) is.element(el, set2)))
  }

  lengthSetsUnique <- sapply(setsUnique, length)
  idxMinSets <- which(lengthSetsUnique == min(lengthSetsUnique))

  remove <- NULL
  for(mS in 1:length(idxMinSets)){
    remove <- c(remove, which(sapply(setsUnique, function(sU) isSubset(setsUnique[idxMinSets[mS]], sU))))
  }

  remove <- setdiff(remove, idxMinSets)
  setsUnique <- if(length(remove) > 0) setsUnique[-remove] else setsUnique

  lengthsSets <- sapply(setsUnique, length)
  setsUnique[order(lengthsSets)]
}


buildTrees <- function(as, tree, addTo){

  occur <- table(unlist(as))
  uniqueVars <- sort(unique(unlist(as)))

  for(i in 1:length(uniqueVars)){
    addTo$AddChild(uniqueVars[i])
  }

  # cur.max <- as.numeric(names(which( occur == max(occur) )))

  # if(length(cur.max) == 1){
  #   addTo$AddChild(cur.max)
  # }else{
  #   for(cm in 1:length(cur.max)){
  #     addTo$AddChild(cur.max[cm])
  #   }
  # }

  tree
}

getListWithIndicesOfDefSets <- function(definingSets, colnamesXX){
  definingSetsColnames <- lapply(definingSets, function(x) colnamesXX[x])
  list(setsUnique = definingSets, definingSetsColnames = definingSetsColnames)
}


getListWithIndicesOfDefSetsOld <- function(definingSets, colnamesXX){
  setsUnique <- NULL

  for(ds in seq_along(definingSets)){
    rootToLeafPaths <- definingSets[[ds]]$Get("pathString", filterFun = function(x) x$isLeaf)
    sets <- lapply(rootToLeafPaths, function(x) sort(as.numeric(strsplit(x, "/")[[1]])))
    setsUnique <- c(setsUnique, sets[which(!duplicated(sets))])
  }

  setsUnique <- setsUnique[which(!duplicated(setsUnique))]
  definingSetsColnames <- lapply(setsUnique, function(x) colnamesXX[x])
  list(setsUnique = setsUnique, definingSetsColnames = definingSetsColnames)
}


