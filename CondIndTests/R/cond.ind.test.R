cond.ind.test <- function(X, Y, Z,
                          method = "KCI",
                          alpha = 0.05,
                          par.test = list(),
                          verbose = FALSE){
  if (method == "KCI"){
    if(!exists('width', where=par.test)) par.test$width <- 0     
    if(!exists('unbiased', where=par.test)) par.test$unbiased <- FALSE     
    if(!exists('approx', where=par.test)) par.test$approx <- TRUE     
    if(!exists('bootstrap', where=par.test)) par.test$bootstrap <- TRUE     
    if(!exists('nRepBs', where=par.test)) par.test$nRepBs <- 500     
    if(!exists('lambda', where=par.test)) par.test$lambda <- 1E-3     
    if(!exists('thresh', where=par.test)) par.test$thresh <- 1E-5     
    if(!exists('numEig', where=par.test)) par.test$numEig <- length(Y)    
    result <- KCI(X, Y, Z,    width = par.test$width,
                  alpha = alpha,
                  unbiased = par.test$unbiased,
                  approx = par.test$approx,
                  bootstrap = par.test$bootstrap,
                  nRepBs = par.test$nRepBs,
                  lambda = par.test$lambda,
                  thresh = par.test$thresh,
                  numEig = par.test$numEig,
                  verbose = verbose)
  } else if (method == "wilcox"){
    result <- NULL 
  } else
  {
    stop('Unknown method.')
  }
  
  # result should be a list containing , e.g. 
  # list(testStatistic = statistic, criticalValue = critVal, pvalue = pVal)
  return(result)
}