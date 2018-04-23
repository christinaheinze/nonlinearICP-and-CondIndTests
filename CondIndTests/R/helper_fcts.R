check_input <- function(X, Y, Z, return_vec = FALSE){
  X <- check_input_single(X, return_vec = return_vec)
  Y <- check_input_single(Y, check_factor = TRUE, return_vec = return_vec)
  Z <- check_input_single(Z, return_vec = return_vec)
}

check_input_single <- function(X, check_factor = FALSE, return_vec = FALSE, str = ""){
  if(is.data.frame(X)){
    if(check_factor) lapply(X, check_E)
    if(return_vec){
      if(NCOL(X) == 1){
        X <- X[,1]
      }else{
        stop(paste("Input", str, "should be one-dimensional"))
      }
    }
  }else if(is.matrix(X)){
    if(check_factor) apply(X, 2, check_E)
    if(return_vec){
      if(NCOL(X) == 1){
        X <- X[,1]
      }else{
        stop(paste("Input", str, "should be one-dimensional"))
      }
    }
  }else if(is.vector(X)){
    if(check_factor) check_E(X)
  }else if(!is.factor(X)){
    stop("check type")
  }
  X
}


check_E <- function(E){
  if(!is.factor(E)){
    uE <- unique(E)
    nruE <- NROW(uE)
    if(nruE < 5){
      warning("E has less than 5 unique values; are you sure that E is not a factor?")
    }
  }
}