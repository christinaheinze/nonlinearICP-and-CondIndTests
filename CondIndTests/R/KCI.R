#' Kernel conditional independence test.
#'
#' @description Tests the null hypothesis that Y and E are independent given X. The distribution of the test
#' statistic under the null hypothesis equals an infinite weighted sum of chi squared variables. This distribution
#' can either be approximated by a gamma distribution or by a Monte Carlo approach. 
#' Hyperparameters can be chosen using Gaussian Process regression.
#'
#' @param Y A vector of length n or a matrix or dataframe with n rows and p columns.
#' @param E A vector of length n or a matrix or dataframe with n rows and p columns.
#' @param X A matrix or dataframe with n rows and p columns.
#' @param width Kernel width; if it is set to zero, the width is chosen automatically (default: 0).
#' @param alpha Significance level (default: 0.05).
#' @param unbiased A boolean variable that indicates whether a bias correction should be applied (default: FALSE).
#' @param gammaApprox A boolean variable that indicates whether the null distribution is approximated by a Gamma
#' distribution. If it is FALSE, a Monte Carlo approach is used (default: TRUE).
#' @param GP Flag whether to use Gaussian Process regression to choose the hyperparameters
#' @param nRepBs Number of draws for the Monte Carlo approach (default: 500).
#' @param lambda Regularization parameter (default: 1e-03).
#' @param thresh Threshold for eigenvalues. Whenever eigenvalues are computed, they are set to zero if they are
#' smaller than thresh times the maximum eigenvalue (default: 1e-05).
#' @param numEig Number of eigenvalues computed (only relevant for computing the distribution under the hypothesis
#' of conditional independence) (default: length(Y)).
#' @param verbose If \code{TRUE}, intermediate output is provided. (default: \code{FALSE}).
#'
#' @return A list with the following entries:
#' \itemize{
#'  \item \code{testStatistic} the statistic Tr(K_{(ddot{(Y)}|X)} * K_{(E|X)})
#'  \item \code{criticalValue} the critical point at the p-value equal to alpha;
#'   obtained by a Monte Carlo approach if \code{gammaApprox = FALSE}, otherwise obtained by Gamma approximation.
#'  \item \code{pvalue} The p-value for the null hypothesis that Y and E are independent given X.
#'  It is obtained by a Monte Carlo approach if \code{gammaApprox = FALSE}, 
#'  otherwise obtained by Gamma approximation.
#'  }
#'
#' @examples
#' # Example 1
#' n <- 100
#' E <- rnorm(n)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * (X)^2 + rnorm(n)
#' KCI(Y, E, X)
#' KCI(Y, X, E)
#'

KCI <- function(Y, E, X,
                width = 0,
                alpha = 0.05,
                unbiased = FALSE,
                gammaApprox = TRUE,
                GP = TRUE,
                nRepBs = 5000,
                lambda = 1E-3,
                thresh = 1E-5,
                numEig = NROW(Y),
                verbose = FALSE){
  
  ## recode input names to match matlab implementation
  Z <- X
  X <- Y
  Y <- E
  
  X <- check_input_single(X)
  dimY <- NCOL(Y)
  if(dimY == 1){
    Y <- check_input_single(Y, return_vec = TRUE, check_factor = TRUE)
  }else{
    Y <- check_input_single(Y, return_vec = FALSE, check_factor = TRUE)
  }
  Z <- check_input_single(Z)
  
  # sample size
  n <- NROW(X)
  
  # normalize the data
  X <- scale(X)
  ZPrime <- scale(Z)
  
  if(any(is.na(ZPrime))){
    Z <- scale(Z, scale = FALSE)
  }else{
    Z <- ZPrime
  }
  rm(ZPrime)
  
  # numbers of variables in Z
  d <- NCOL(Z)
  
  # kernel width
  if(width == 0){
    if(n <= 200){
      width <- 1.2 # 0.8
    }else if(n < 1200){
      width <- 0.7 # 0.5
    }else{
      width <- 0.4 # 0.3
    }
  }
  
  kernPrecision <- 1/(width^2 * d) # width = sqrt(2) sigma
  
  # for centering of the data in feature space
  H <- diag(n) - matrix(1,n,n)/n
  
  # compute kernels for (X, Z) (which is called Kx in the Matlab implementation)
  KXZ <- rbfKernel1(cbind(X, Z/2),c (kernPrecision,1))$kx
  # centralized kernel matrix
  # KYZ <- H %*% KYZ %*% H
  KXZ <- crossprod(H, KXZ) %*% H
  
  if(is.factor(Y) & dimY == 1){
    # delta kernel for categorical variable E
    Ynum <- as.numeric(Y)
    KY <- sapply(Ynum, function(i) i == Ynum) %*% diag(n)
  }else{
    if(is.data.frame(Y) & all(sapply(Y, is.factor))){
      Y <- as.matrix(data.frame(lapply(Y, function(i) as.numeric(as.character(i)))))
    }
    Y <- scale(Y)
    KY <- rbfKernel1(Y, c(kernPrecision,1))$kx
  }
  # centralized kernel matrix
  KY <- crossprod(H, KY) %*% H
  
  if (GP) {
    KX <- KXZ
    
    # compute eigen vals and vecs
    numEigKX <- min(400, floor(n/4))
    eigenKX <- svd((KX + t(KX))/2, nu = numEigKX, nv = 0)
    eigValsKX <- eigenKX$d[1:numEigKX]
    eigVecsKX <- eigenKX$u
    
    numEigKY <- min(200, floor(n/5))
    eigenKY <- svd((KY+t(KY))/2, nu = numEigKY, nv = 0)
    eigValsKY <- eigenKY$d[1:numEigKY]
    eigVecsKY <- eigenKY$u
    
    # filter using condition:
    index_KX <- which(eigValsKX > max(eigValsKX) * thresh)
    index_KY <- which(eigValsKY > max(eigValsKY) * thresh)
    eigValsKX <- eigValsKX[index_KX]
    eigValsKY <- eigValsKY[index_KY]
    eigVecsKX <- eigVecsKX[,index_KX]
    eigVecsKY <- eigVecsKY[,index_KY]
    
    # Find optimal hyperparamerters
    
    # initial hyperparam values
    logtheta0 <- c(log(width * sqrt(d)) * c(rep(1, d)), 0, log(sqrt(0.1)))
    
    logtheta_x <- logtheta0
    logtheta_y <- logtheta0
    
    # ifelse statements due to automatic casting when
    # eigVals* has length 1
    if(length(eigValsKX) == 1){
      gpr_KX <- 2 * sqrt(n) * matrix(eigVecsKX) %*% matrix(sqrt(eigValsKX) / sqrt(eigValsKX[1]))
    }else{
      gpr_KX <- 2 * sqrt(n) * eigVecsKX %*% diag(sqrt(eigValsKX)) / sqrt(eigValsKX[1])
    }
    if(length(eigValsKY) == 1){
      gpr_KY <- 2 * sqrt(n) * matrix(eigVecsKY) %*% matrix(sqrt(eigValsKY) / sqrt(eigValsKY[1]))
    }else{
      gpr_KY <- 2 * sqrt(n) * eigVecsKY %*% diag(sqrt(eigValsKY)) / sqrt(eigValsKY[1])
    }
    
    covfunc <- cov.SEard_with_noise
    
    # get optimal hyperparamerters using GP regression
    tryCatch({
      logtheta_x <- gpr.find_optimal_hyperparams(logtheta0, Z, gpr_KX, covfunc, verbose)
      logtheta_y <- gpr.find_optimal_hyperparams(logtheta0, Z, gpr_KY, covfunc, verbose)
    }, error = function(e) {
      warning(
        "Finding optimal kernel hyperparams failed with error: \"",
        e[1],
        "\" Continuing with heuristic hyperparams."
      )
    })
    
    # Note: in the conditional case, no need to do centering, as the regression
    # will automatically enforce that.
    
    # Build the KZ kernel matrixes with correct hyperparams
    KZ_x <- cov.SEard$gram(head(logtheta_x, n=d+1), Z)
    KZ_y <- cov.SEard$gram(head(logtheta_y, n=d+1), Z)
    
    # Kernel matrices of the errors
    RZ_x <- diag(n) - KZ_x %*% MASS::ginv(KZ_x + exp(2*tail(logtheta_x, n=1))*diag(n)) # Matlab: P1_x
    RZ_y <- diag(n) - KZ_y %*% MASS::ginv(KZ_y + exp(2*tail(logtheta_y, n=1))*diag(n)) # Matlab: P1_y
    
    KXZ <- RZ_x %*% tcrossprod(KX, RZ_x)
    KYZ <- RZ_y %*% tcrossprod(KY, RZ_y)
    
    statistic <- sum(diag(KXZ %*% KYZ))
    
    # degrees of freedom
    dfY <- sum(diag(diag(n) - RZ_y))
    dfX <- sum(diag(diag(n) - RZ_x))
    
  } else {
    # kernel for conditioning set Z
    KZ <- rbfKernel1(Z, c(kernPrecision,1))$kx
    # centralized kernel matrix
    # KZ <- H %*% KZ %*% H
    KZ <- crossprod(H, KZ) %*% H
    
    # kernel matrices of the errors
    RZ <- diag(n) - KZ%*%MASS::ginv(KZ + lambda*diag(n)) #expensive
    # RZ <- diag(n) - KZ%*%pseudoinverse(KZ + lambda*diag(n)) #expensive
    # KYZ <- RZ %*% KYZ %*% t(RZ) #Eq. (11)
    KXZ <- RZ %*% tcrossprod(KXZ,RZ) #Eq. (11)
    
    # KEZ <- RZ %*% KE %*% t(RZ) # Eq. (12)
    KYZ <- RZ %*% tcrossprod(KY,RZ) # Eq. (12)
    
    # calculate the statistic
    statistic <- sum(diag(KXZ %*% KYZ))
    
    # degrees of freedom
    dfY <- dfX <- sum(diag(diag(n) - RZ))
  }
  
  # calculate the eigenvalues
  # due to numerical issues, KXZ and KYZ may not be symmetric:
  eigenXZ <- svd((KXZ+t(KXZ))/2, nu = numEig, nv = 0)
  eigValsKXZ <- eigenXZ$d[1:numEig]
  eigVecsKXZ <- eigenXZ$u
  eigenYZ <- svd((KYZ+t(KYZ))/2, nu = numEig, nv = 0)
  eigValsKYZ <- eigenYZ$d[1:numEig]
  eigVecsKYZ <- eigenYZ$u
  
  # calculate the product of the square root of the eigvector and the eigenvector
  IIxz <- which(eigValsKXZ > max(eigValsKXZ) * thresh)
  IIyz <- which(eigValsKYZ > max(eigValsKYZ) * thresh)
  eigValsKXZ <- eigValsKXZ[IIxz]
  eigVecsKXZ <- eigVecsKXZ[,IIxz]
  eigValsKYZ <- eigValsKYZ[IIyz]
  eigVecsKYZ <- eigVecsKYZ[,IIyz]
  
  if(length(eigValsKXZ) == 1){
    eivProdXZ <- matrix(eigVecsKXZ * sqrt(eigValsKXZ), ncol = 1)
  }else{
    eivProdXZ <- tcrossprod(eigVecsKXZ, diag(sqrt(eigValsKXZ)))
  }
  
  if(length(eigValsKYZ) == 1){
    eivProdYZ <- matrix(eigVecsKYZ * sqrt(eigValsKYZ), ncol = 1)
  }else{
    eivProdYZ <- tcrossprod(eigVecsKYZ, diag(sqrt(eigValsKYZ)))
  }
  
  rm(eigValsKXZ, eigVecsKXZ, eigValsKYZ, eigVecsKYZ)
  
  # calculate their product
  nEivProdXZ <- ncol(eivProdXZ)
  nEivProdYZ <- ncol(eivProdYZ)
  
  sizeU <- nEivProdXZ*nEivProdYZ
  uu <- matrix(0, n, sizeU)
  
  for(i in 1:nEivProdXZ){
    for(j in 1:nEivProdYZ){
      uu[,((i-1)*nEivProdYZ + j)] <- eivProdXZ[,i] * eivProdYZ[,j]
    }
  }
  
  if(sizeU > n){
    # uuProd <- uu %*% t(uu)
    uuProd <- tcrossprod(uu)
    
  }else{
    # uuProd <- t(uu) %*% uu
    uuProd <- crossprod(uu)
  }
  if(!gammaApprox){ #if bootstrap
    eigUU <- eigen(uuProd, only.values = TRUE)$values
    keepN <- min(n,sizeU)
    eigUU <- eigUU[1:keepN]
    
    IIf <- which(eigUU > max(eigUU) * thresh)
    eigUU <- eigUU[IIf]
  }
  critVal <- NULL
  pVal <- NULL
  
  if(!gammaApprox){ #if bootstrap
    if(length(eigUU) * n < 1E6){
      fRand1 <- matrix(rchisq(n = length(eigUU) * nRepBs, df = 1),
                       nrow = length(eigUU), ncol = nRepBs)
      
      if(unbiased){
        nullDistr <- n^2/(n-1-dfX)/(n-1-dfY) * crossprod(eigUU, fRand1)
      }else{
        nullDistr <- crossprod(eigUU, fRand1)
      }
      
    }else{
      # iteratively calcuate the null dstr to save memory
      nullDistr <- matrix(0, 1, nRepBs)
      len <- max(floor(1E6/n), 100)
      maxIter <- floor(length(eigUU)/len)
      
      for(iter in 1:maxIter){
        fRand1 <- matrix(rchisq(n = len * nRepBs, df = 1),
                         nrow = len, ncol = nRepBs)
        if(unbiased){
          nullDistr <- nullDistr +
            n^2/(n-1-dfX)/(n-1-dfY) * crossprod(eigUU[((iter-1)*len+1):(iter*len)], fRand1)
        }else{
          nullDistr <- nullDistr + crossprod(eigUU[((iter-1)*len+1):(iter*len)], fRand1)
        }
      }
    }
    sortNullDistr <- sort(nullDistr)
    critVal <- sortNullDistr[ceiling((1-alpha)*nRepBs)]
    pVal <- sum(nullDistr > statistic)/nRepBs
  }else{
    
    meanApprox <- sum(diag(uuProd))
    varApprox <- 2*sum(diag(uuProd%*%uuProd))
    kApprox <- meanApprox^2/varApprox
    kernPrecisionApprox <- varApprox/meanApprox
    critVal <- qgamma(1-alpha,
                      shape = kApprox,
                      scale = kernPrecisionApprox,
                      lower.tail = TRUE)
    pVal <- 1 - pgamma(statistic,
                       shape = kApprox,
                       scale = kernPrecisionApprox,
                       lower.tail = TRUE)
  }
  
  list(testStatistic = statistic, criticalValue = critVal, pvalue = pVal)
}
