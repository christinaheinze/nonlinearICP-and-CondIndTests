rbfKernel <- function(x, xKern, theta){

  precision <- theta[1]
  scaling <- theta[2]

  distMat <- dist2(x, xKern)

  if(precision == 0){
    precision <- 2/median(distMat[lower.tri(distMat) > 0])
  }

  precisionDivBy2 <- precision/2
  kx <- scaling*exp(-precisionDivBy2*distMat)
  bwNew <- 1/precision

  list(kx = kx, bwNew = bwNew)
}

rbfKernel1 <- function(x, theta){

  precision <- theta[1]
  scaling <- theta[2]

  distMat <- dist1(x)

  if(precision == 0){
    precision <- 2/median(distMat[lower.tri(distMat) > 0])
  }

  precisionDivBy2 <- precision/2
  kx <- scaling*exp(-precisionDivBy2*distMat)
  bwNew <- 1/precision

  list(kx = kx, bwNew = bwNew)
}