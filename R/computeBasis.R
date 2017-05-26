computeBasis <- function(basis, X, p, degree, nSub, verbose){
  if(basis == "polynomial"){
    XBasis <- polyBasis(X, p, degree)

    fctBasisExpansion <- function(Xprime){
      polyBasis(Xprime, p, degree)
    }

    if(verbose){
      cat(paste("\nPolynomial basis has been created. Design matrix has",
                NCOL(XBasis), "columns. Degree", degree))
    }

  }else if(is.element(basis, c("nystrom", "nystrom_poly", "fourier"))){


    fctBasisExpansion <- randomFeatsFct(as.matrix(X), nSub, type = basis, degree = degree)
    XBasis <- fctBasisExpansion(as.matrix(X))

    if(verbose){
      cat(paste("\nRandom features with", basis, "basis have been created.
                Design matrix has",
                NCOL(XBasis), "columns. Degree", degree))
    }

  }else{
    stop("\nEstimation procedure not supported.")
  }

  list(XBasis = XBasis, fctBasisExpansion = fctBasisExpansion)
}


randomFeatsFct <- function(x,k,type="fourier", degree = 2) {
  s <- sigest(x,scaled=NULL)[2]

  if(type == "linear") {
    return(function(x0) x0)
  }

  if(type == "nystrom") {
    w <- x[sample(1:NROW(x),k),]
    return(function(x0) kernelMatrix(rbfdot(s),x0,w))
  }

  if(type == "nystrom_poly") {
    w <- x[sample(1:NROW(x),k), , drop = FALSE]
    return(function(x0) kernelMatrix(polydot(degree = degree),x0,w))
  }

  if(type == "fourier") {
    w <- matrix(rnorm(NCOL(x)*k,sd=sqrt(2*s)),NCOL(x))
    b <- runif(k,0,2*pi)
    f <- function(x0) x0%*%w+t(matrix(b,k,NROW(x0)))
    return(function(x0) cos(f(x0)))
  }
}

polyBasis <- function(X, p, q){

  if(is.null(colnames(X))){
    X <- as.matrix(X)
    colnames(X) <- 1:p
  }

  XBasis <- matrix(nrow=NROW(X),ncol=p*q)
  colX <- character(2*p)
  cc <- 0
  for (k in 1:p){
    cc <- cc+1
    XBasis[,cc] <- X[,k]
    colX[cc] <- colnames(X)[k]
    for (qc in 2:q){
      cc <- cc+1
      XBasis[,cc] <- X[,k]^qc
      colX[cc] <- paste(colnames(X)[k],"_",qc,sep="")
    }
  }
  colnames(XBasis) <- colX
  XBasis
}
