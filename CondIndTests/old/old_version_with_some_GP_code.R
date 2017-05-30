#   To test if Y and E (discrete) are independent, given X.
#   INPUT:
#   The number of rows of x and y is the sample size.
#   alpha is the significance level (we suggest 1%).
#   width contains the kernel width.
#   Output:
#   critVal: the critical point at the p-value equal to alpha obtained by bootstrapping.
#   statistic: the statistic Tr(K_{\ddot{X}|Z} * K_{Y|Z}).
#   pVal: the p value obtained by bootstrapping.
#   critValAppr: the critical value obtained by Gamma approximation.
#   p_apppr: the p-value obtained by Gamma approximation.
#   If Sta > critVal, the null hypothesis (x is independent from y) is rejected.
#


KCIDiscrete <- function(Y, E, X, width = 0, alpha = 0.01, unbiased = FALSE,
                        GP = FALSE, approx = TRUE, bootstrap = TRUE,
                        nRepBs = 5000, lambda = 1E-3, thresh = 1E-5,
                        numEig = length(Y), nIterGP = 350, verbose = FALSE){

  # sample size
  n <- length(Y)

  # normalize the data
  Y <- scale(Y)
  E <- scale(E) #TODO: correct to scale E?
  X <- scale(X)
  # numbers of variables in X
  d <- ncol(X)

  # find number of samples from each environment (atm: just two envs supported!)
  nInE <- as.data.frame(table(E))
  n1 <- nInE[1,"Freq"]
  n2 <- nInE[2,"Freq"]

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
  # KYX <- rbfKernel1(cbind(Y, X/2),c(kernPrecision,1))$kx
  KYX <- rbfKernel(cbind(Y, X/2),cbind(Y, X/2),c(kernPrecision,1))$kx
  
  
  # centralized kernel matrix
  KYX <- H %*% KYX %*% H

  # delta kernel for discrete variable E
  KE <- rbfKernel(E, E, c(kernPrecision,1))$kx #labelDot(n1,n2) #TODO: which one to use??
  # KE <- rbfKernel1(E, c(kernPrecision,1))$kx 
  KE <- H %*% KE %*% H
  # KE <- discrete_grammat(E)

  if(GP){
    stop("Implementation using Gaussian processes incomplete!")
    ### learning the hyperparameters
    if(verbose) cat('\nOptimizing hyperparameters in GP regression...\n')

    options <- gpOptions()
    options$kern$comp = list("rbf", "white")
    YX <- cbind(Y, X)
    modelXY <- gpCreate(ncol(X), ncol(YX), X, YX, options)
    optimOutXY <- gpOptimise(modelXY, ifelse(verbose, 1, 0), nIterGP, FALSE)
    parsXY <- optimOutXY$params$xmin
#     inverseWidth <- pars[1]
#     signalVariance <- pars[2]
    noiseVarianceXY <- parsXY[3]

    modelXE <- gpCreate(ncol(X), 1, X, matrix(E, length(E), 1), options)
    optimOutXE <- gpOptimise(modelXE, ifelse(verbose, 1, 0), nIterGP, FALSE)
    parsXE <- optimOutXE$params$xmin
    noiseVarianceXE <- parsXE[3]

    if(verbose) cat(paste("\nNoise variance XY", round(exp(noiseVarianceXY), 4),
                          "\nNoise variance XE", round(exp(noiseVarianceXE), 4)))

#     optionsX = gpOptions()
#     optionsX$kern$comp = list("rbf")
#     modelX <- gpCreate(ncol(X), ncol(X), X, X, optionsX)
#
#     kernXY <- gpExpandParam(modelX, parsXY)$K_uu
#     kernXE <- gpExpandParam(modelX, parsXE)$K_uu

    kernXY <- rbfKernel(X, X, exp(parsXY[1:2]))$kx
    kernXE <- rbfKernel(X, X, exp(parsXE[1:2]))$kx

    # Kernel matrices of the errors
    RY <-  diag(n) - kernXY%*%MASS::ginv(kernXY + exp(noiseVarianceXY)*diag(n))
    KYX <- RY %*% KYX %*% t(RY)
    RE <- diag(n) - kernXE%*%MASS::ginv(kernXE + exp(noiseVarianceXE)*diag(n))
    KEX <- RE %*% KE %*% t(RE)
    # Jonas: das verstehe ich nicht ganz. sollten die R_Z's (Kun) = R_X's (hier) in (11) und (12) nicht identisch sein? liegt das am GP? naja, unten ist es ja anders...
    #

    # calculate the statistic
    statistic <- sum(diag(KYX %*% KEX))/n

    # degrees of freedom
    dfY <- sum(diag(diag(n) - RY))
    dfE <- sum(diag(diag(n) - RE))

  }else{
    # kernel for conditioning set X
    # KX <- rbfKernel1(X, c(kernPrecision,1))$kx
    KX <- rbfKernel(X, X,c(kernPrecision,1))$kx
    
    # centralized kernel matrix
    KX <- H %*% KX %*% H

    # kernel matrices of the errors
    RX <- diag(n) - KX%*%MASS::ginv(KX + lambda*diag(n))
    KYX <- RX %*% KYX %*% t(RX) #Eq. (11)
    KEX <- RX %*% KE %*% t(RX) # Eq. (12)
    # Jonas: es ist zugegebener Weise etwas komisch bei der Regression von E auf Z(Kun)=X(hier) eine ridge regression zu machen, aber ich wuerde einfach einen Kommentar einfuegen und fertig :-).
    #

    # calculate the statistic
    statistic <- sum(diag(KYX %*% KEX))

    # degrees of freedom
    dfE <- dfY <- sum(diag(diag(n) - RX))

  }

  # calculate the eigenvalues
  # due to numerical issues, KYX and KEX may not be symmetric:
  eigenYX <- svd((KYX+t(KYX))/2, nu = numEig, nv = 0) #eigen((KYX+t(KYX))/2)
  eigValsKYX <- eigenYX$d[1:numEig] #eigenYX[1]
  eigVecsKYX <- eigenYX$u #eigenYX[2]
  eigenEX <- svd((KEX+t(KEX))/2, nu = numEig, nv = 0) #eigen((KEX+t(KEX))/2)
  eigValsKEX <- eigenEX$d[1:numEig] #eigenEX[1]
  eigVecsKEX <- eigenEX$u #eigenEX[2]


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
    eivProdYX <- eigVecsKYX %*% diag(sqrt(eigValsKYX))
  }

  if(length(eigValsKEX) == 1){
    eivProdEX <- matrix(eigVecsKEX * sqrt(eigValsKEX), ncol = 1)
  }else{
    eivProdEX <- eigVecsKEX %*% diag(sqrt(eigValsKEX))
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
    uuProd <- uu %*% t(uu)
  }else{
    uuProd <- t(uu) %*% uu
  }
  if(bootstrap){
    eigUU <- eigen(uuProd, only.values = TRUE)$values
    keepN <- min(n,sizeU)
    eigUU <- eigUU[1:keepN]

    IIf <- which(eigUU > max(eigUU) * thresh)
    eigUU <- eigUU[IIf]
  }
  critVal <- NULL
  pVal <- NULL

  if(bootstrap){
    if(length(eigUU) * n < 1E6){
      fRand1 <- matrix(rchisq(n = length(eigUU) * nRepBs, df = 1),
                       nrow = length(eigUU), ncol = nRepBs)

      if(unbiased){
        nullDistr <- n^2/(n-1-dfY)/(n-1-dfE) * t(eigUU) %*% fRand1 ### TODO: dfs / ??
      }else{
        nullDistr <- t(eigUU) %*% fRand1
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
            n^2/(n-1-dfY)/(n-1-dfE) * t(eigUU[((iter-1)*len+1):(iter*len)]) %*% fRand1 # TODO: dfs / ??
        }else{
          nullDistr <- nullDistr + t(eigUU[((iter-1)*len+1):(iter*len)]) %*% fRand1
        }
      }
    }
    #nullDistr <- nullDistr/n
    sortNullDistr <- sort(nullDistr)
    critVal <- sortNullDistr[ceiling((1-alpha)*nRepBs)]
    pVal <- sum(nullDistr > statistic)/nRepBs
  }

  critValAppr <- NULL
  pValAppr <- NULL

  if(approx){
    meanApprox <- sum(diag(uuProd))
    varApprox <- 2*sum(diag(uuProd^2))
    kApprox <- meanApprox^2/varApprox
    kernPrecisionApprox <- varApprox/meanApprox
    critValAppr <- qgamma(1-alpha,
                          shape = kApprox,
                          scale = kernPrecisionApprox,
                          lower.tail = TRUE)
    pValAppr <- 1 - pgamma(statistic,
                           shape = kApprox,
                           scale = kernPrecisionApprox,
                           lower.tail = TRUE)
  }

  list(testStatistic = statistic, criticalValue = critVal,
       pValue = pVal, criticalValueApprox = critValAppr, pValueApprox = pValAppr)
}
