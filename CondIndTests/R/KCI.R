#' Kernel conditional independence test.
#'
#' @description Tests the null hypothesis that Y and E are independent given X. The distribution of the test
#' statistic under the null hypothesis equals an infinite weighted sum of chi squared variables. This distribution
#' can either be approximated by a gamma distribution or by a Monte Carlo approach. Choosing the hyperparameters
#' by Gaussian Process regression is not yet implemented.
#'
#' @param Y A vector of length n or a matrix or dataframe with n rows and p columns.
#' @param E A vector of length n or a matrix or dataframe with n rows and p columns.
#' @param X A matrix or dataframe with n rows and p columns.
#' @param width Kernel width; if it is set to zero, the width is chosen automatically (default: 0).
#' @param alpha Significance level (default: 0.05).
#' @param unbiased A boolean variable that indicates whether a bias correction should be applied (default: FALSE).
#' @param gammaApprox A boolean variable that indicates whether the null distribution is approximated by a Gamma
#' distribution. If it is FALSE, a Monte Carlo approach is used (default: TRUE).
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
#'  \item \code{testStatistic} the statistic Tr(K_{(ddot{(X)}|Z)} * K_{(Y|Z)})
#'  \item \code{criticalValue} the critical point at the p-value equal to alpha;
#'   obtained by a Monte Carlo approach if \code{gammaApprox = FALSE}, otherwise obtained by Gamma approximation.
#'  \item \code{pvalue} The p-value for the null hypothesis that Y and E are independent given X.
#'  It is obtained by a Monte Carlo approach if \code{gammaApprox = FALSE}, otherwise obtained by Gamma approximation.
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
                nRepBs = 500,
                lambda = 1E-3,
                thresh = 1E-5,
                numEig = if(is.data.frame(Y) | is.matrix(Y)) nrow(Y) else if(is.vector(Y)) length(Y),
                verbose = FALSE){

  if(!is.factor(E)){
    uE <- unique(E)
    nruE <- if(is.data.frame(E) | is.matrix(E)) nrow(uE) else if(is.vector(E)) length(uE)
    if(nruE < 5){
      warning("E has less than 5 unique values; are you sure that E is not a factor?")
    }
  }

  # sample size
  n <- if(is.data.frame(Y) | is.matrix(Y)) nrow(Y) else if(is.vector(Y)) length(Y)

  # normalize the data
  Y <- scale(Y)
  XPrime <- scale(X)

  if(any(is.na(XPrime))){
    X <- scale(X, scale = FALSE)
  }else{
    X <- XPrime
  }
  rm(XPrime)

  # numbers of variables in X
  d <- ncol(X)

  # kernel width
  if(width == 0){
    if(n <= 200){
      width <- 0.8
    }else if(n < 1200){
      width <- 0.5
    }else{
      width <- 0.3
    }
  }

  kernPrecision <- 1/(width^2 * d) # width = sqrt(2) sigma

  # for centering of the data in feature space
  H <- diag(n) - matrix(1,n,n)/n
  # compute kernels for (Y, X)
  KYX <- rbfKernel1(cbind(Y, X/2),c(kernPrecision,1))$kx
  # centralized kernel matrix
  # KYX <- H %*% KYX %*% H
  KYX <- crossprod(H, KYX) %*% H

  if(is.factor(E)){
    # show(E)
    # delta kernel for discrete variable E
    Enum <- as.numeric(E)
    KE <- (Enum^2 == (Enum %*% t(Enum))) %*% diag(n)
  } else {
    E <- scale(E)
    KE <- rbfKernel1(E, c(kernPrecision,1))$kx
    # centralized kernel matrix
    # KE <- H %*% KE %*% H
  }
  KE <- crossprod(H, KE) %*% H

  # kernel for conditioning set X
  KX <- rbfKernel1(X, c(kernPrecision,1))$kx
  # centralized kernel matrix
  # KX <- H %*% KX %*% H
  KX <- crossprod(H, KX) %*% H

  # kernel matrices of the errors
  RX <- diag(n) - KX%*%MASS::ginv(KX + lambda*diag(n)) #expensive
  # RX <- diag(n) - KX%*%pseudoinverse(KX + lambda*diag(n)) #expensive
  # KYX <- RX %*% KYX %*% t(RX) #Eq. (11)
  KYX <- RX %*% tcrossprod(KYX,RX) #Eq. (11)

  # KEX <- RX %*% KE %*% t(RX) # Eq. (12)
  KEX <- RX %*% tcrossprod(KE,RX) # Eq. (12)

  # calculate the statistic
  statistic <- sum(KYX * t(KEX)) #sum(diag(KYX %*% KEX))

  # degrees of freedom
  dfE <- dfY <- sum(diag(diag(n) - RX))

  # calculate the eigenvalues
  # due to numerical issues, KYX and KEX may not be symmetric:
  eigenYX <- svd((KYX+t(KYX))/2, nu = numEig, nv = 0)
  eigValsKYX <- eigenYX$d[1:numEig]
  eigVecsKYX <- eigenYX$u
  eigenEX <- svd((KEX+t(KEX))/2, nu = numEig, nv = 0)
  eigValsKEX <- eigenEX$d[1:numEig]
  eigVecsKEX <- eigenEX$u

  # calculate the product of the square root of the eigvector and the eigenvector
  IIyx <- which(eigValsKYX > max(eigValsKYX) * thresh)
  IIex <- which(eigValsKEX > max(eigValsKEX) * thresh)
  eigValsKYX <- eigValsKYX[IIyx]
  eigVecsKYX <- eigVecsKYX[,IIyx]
  eigValsKEX <- eigValsKEX[IIex]
  eigVecsKEX <- eigVecsKEX[,IIex]

  if(length(eigValsKYX) == 1){
    eivProdYX <- matrix(eigVecsKYX * sqrt(eigValsKYX), ncol = 1)
  }else{
    eivProdYX <- tcrossprod(eigVecsKYX, diag(sqrt(eigValsKYX)))
  }

  if(length(eigValsKEX) == 1){
    eivProdEX <- matrix(eigVecsKEX * sqrt(eigValsKEX), ncol = 1)
  }else{
    eivProdEX <- tcrossprod(eigVecsKEX, diag(sqrt(eigValsKEX)))
  }

  rm(eigValsKYX, eigVecsKYX, eigValsKEX, eigVecsKEX)

  # calculate their product
  nEivProdYX <- ncol(eivProdYX)
  nEivProdEX <- ncol(eivProdEX)

  sizeU <- nEivProdYX*nEivProdEX
  uu <- matrix(0, n, sizeU)

  for(i in 1:nEivProdYX){
    for(j in 1:nEivProdEX){
      uu[,((i-1)*nEivProdEX + j)] <- eivProdYX[,i] * eivProdEX[,j]
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
        nullDistr <- n^2/(n-1-dfY)/(n-1-dfE) * crossprod(eigUU, fRand1)
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
            n^2/(n-1-dfY)/(n-1-dfE) * crossprod(eigUU[((iter-1)*len+1):(iter*len)], fRand1)
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
    varApprox <- 2*sum(diag(uuProd^2))
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
