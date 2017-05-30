# %DIST2	Calculates squared distance between two sets of points.
# %
# %	Description
# %	D = DIST2(X, C) takes two matrices of vectors and calculates the
# %	squared Euclidean distance between them.  Both matrices must be of
# %	the same column dimension.  If X has M rows and N columns, and C has
# %	L rows and N columns, then the result has M rows and L columns.  The
# %	I, Jth entry is the  squared distance from the Ith row of X to the
# %	Jth row of C.
# %
#
# %	Copyright (c) Ian T Nabney (1996-2001)
# % All rights reserved.  See the file COPYING for license terms.

dist2 <- function(X, C){
  dimsX <- dim(X)
  nX <- dimsX[1]
  pX <- dimsX[2]
  dimsC <- dim(C)
  nC <- dimsC[1]
  pC <- dimsC[2]

  if(pX != pC){
    stop('Data dimension does not match dimension of centres')
  }


  sqEucDis <- function(a,b){
    sum((a-b)^2)
  }

  sqEucDisVecMat <- function(vec, mat){
    apply(mat, 1, function(i) sqEucDis(i, vec))
  }

  distMat <- t(apply(X, 1, function(i) sqEucDisVecMat(i, C)))

  # Rounding errors occasionally cause negative entries in distMat
  if(any(distMat < 0)){
    distMat[distMat < 0] <- 0
  }

  distMat
}

dist1 <- function(X){
  xnorm <- as.matrix(dist(X,method="euclidean",diag=TRUE,upper=TRUE))
  xnorm^2
}
