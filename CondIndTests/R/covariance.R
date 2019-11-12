#' R implemantation of covariance functions from vers. 2 of the GPML toolbox by
#' Carl Edward Rasmussen and Christopher K. I. Williams distributed with their
#' book "Gaussian Processes for Machine Learning", the MIT Press 2006.
#'
#' Copyright (c) 2005-2018 Carl Edward Rasmussen & Hannes Nickisch. All rights
#' reserved.
#'
#' Original code was distibuted under the FreeBSD license. See
#' http://www.gaussianprocess.org/gpml/code/matlab/Copyright for additional
#' information.
#'
#' Specificly, this file implements he squared exponentail covariance function
#' with both isotrophic (iso) and automatic relevance determination (ard)
#' distance measures with and without "white noise".
#'

#' Squared Exponential covariance function with automatic relevance
#' determination (ard) distance measure
cov.SEard <- list(
  gram = function(logtheta, x) { cov.gram.SEard(logtheta, x) },
  deriv = function(logtheta, x, z) { cov.deriv.SEard(logtheta, x, z) },
  # one length scale per dimension and a signal variance
  numParams = function(x) { dim(x)[2] + 1 }
)

#' Squared Exponential covariance function with isotrophic (iso) distance
#' measure
cov.SEiso <- list(
  gram = function(logtheta, x) { cov.gram.SEiso(logtheta, x) },
  deriv = function(logtheta, x, z) { cov.deriv.SEiso(logtheta, x, z) },
  # one length scale and signal variance
  numParams = function(x) { 2 }
)

#' Squared Exponential covariance function with isotrophic (iso) distance
#' measure with "white noise"
cov.SEard_with_noise <- list(
  gram = function(logtheta, x) {
    cov.gram.sum(cov.SEard, cov.noise, logtheta, x)
  },
  deriv = function(logtheta, x, z) {
    cov.deriv.sum(cov.SEard, cov.noise, logtheta, x, z)
  },
  # one length scale per dimension, signal variance and a noise variance
  numParams = function(x) { dim(x)[2] + 2 }
)

#' Squared Exponential covariance function with isotrophic (iso) distance
#' measure with "white noise".
cov.SEiso_with_noise <- list(
  gram = function(logtheta, x) {
    cov.gram.sum(cov.SEiso, cov.noise, logtheta, x)
  },
  deriv = function(logtheta, x, z) {
    cov.deriv.sum(cov.SEiso, cov.noise, logtheta, x, z)
  },
  # one length scale, signal variance and a noise variance
  numParams = function(x) { 3 }
)

#' "White noise" covariance functions
cov.noise <- list(
  gram = function(logtheta, x) { cov.gram.noise(logtheta, x) },
  deriv = function(logtheta, x, z) { cov.deriv.noise(logtheta, x, z) },
  numParams = function(x) { 1 }
)

#' Compute the covariance matrix of Squared Exponential covariance function
#' with isotropic (iso) distance measure
#' @param logtheta vector of log hyperparameters.
#' @param x matrix of training inputs of size n x D.
#' @return covariance matrix.
cov.gram.SEiso <- function(logtheta, x) {
  ell <- exp(logtheta[1])       # characteristic length scale
  sf2 <- exp(2 * logtheta[2])   # signal variance
  m <- t(x) * (1 / ell)

  sf2 * exp(-sq_dist(m,m) / 2)
}

#' Compute the derivative matrix of Squared Exponential covariance function
#' with isotropic (iso) distance measure
#' @param logtheta vector of log hyperparameters.
#' @param x matrix of training inputs of size n x D.
#' @param z index of hyperparam in logtheta to compute deriv. with respect to.
#' @return derivative matrix.
cov.deriv.SEiso <- function(logtheta, x, z) {
  if (length(z) > 1) {
    stop("z must be a number type")
  }

  ell <- exp(logtheta[1])       # characteristic length scale
  sf2 <- exp(2 * logtheta[2])   # signal variance
  m <- t(x) * (1 / ell)

  if (z == 1) { # first parameter
    A <- sf2 * exp(-sq_dist(m,m) / 2) * sq_dist(m, m)
  } else {      # second parameter
    A <- 2 * sf2 * exp(-sq_dist(m,m) / 2)
  }

  A
}

#' Compute the covariance matrix of Squared Exponential covariance function
#' with automatic relevance determination (ard) distance measure.
#' @param logtheta vector of log hyperparameters.
#' @param x matrix of training inputs of size n x D.
#' @return covariance matrix.
cov.gram.SEard <- function(logtheta, x) {
  D <- dim(x)[2]                # number of ARD parameters
  ell <- exp(logtheta[1:D])     # characteristic length scale
  sf2 <- exp(2*logtheta[D+1])   # signal variance

  if (is.infinite(sf2)) {
    stop("`sf2` must be finite. Signal variance is probably to large.")
  }

  m <- diag(1/ell, D) %*% t(x)
  sf2 * exp(-sq_dist(m,m) / 2)
}

#' Compute the derivative matrix of Squared Exponential covariance function
#' with automatic relevance determination (ard) distance measure.
#' @param logtheta vector of log hyperparameters.
#' @param x matrix of training inputs of size n x D.
#' @param z index of hyperparam in logtheta to compute deriv. with respect to.
#' @return derivative matrix.
cov.deriv.SEard <- function(logtheta, x, z) {
  if (length(z) > 1) {
    stop("z must be a number type")
  }

  D <- dim(x)[2]                  # number of ARD parameters
  ell <- exp(logtheta[1:D])       # characteristic length scale
  sf2 <- exp(2*logtheta[D + 1])   # signal variance

  m <- diag(1/ell, D) %*% t(x)
  K <- sf2 * exp(-sq_dist(m, m) / 2)
  if (z <= D) {                   # length scale parameters
    m_sub <- t(x[,1:z]) * (1 / ell)
    A <- K * sq_dist(m_sub, m_sub)
  } else {                        # magnitude parameter
    A <- 2 * K
  }

  A
}

#' Compute the covariance matrix of "white noise" (covariance) function to be
#' used with eg. `cov.sum`.
#' @param logtheta vector of log hyperparameters.
#' @param x matrix of training inputs of size n x D.
#' @return covariance matrix.
cov.gram.noise <- function(logtheta = NULL, x) {
  s2 <- exp(2 * as.matrix(logtheta)[1])          # noise variance
  res = s2 * diag(dim(x)[1])
}

#' Compute the derivative matrix of "white noise" (covariance) function to be
#' used with eg. `cov.sum`.
#' @param logtheta vector of log hyperparameters.
#' @param x matrix of training inputs of size n x D.
#' @param z index of hyperparam in logtheta to compute deriv. with respect to.
#' @return derivative matrix.
cov.deriv.noise <- function(logtheta = NULL, x, z) {
  s2 <- exp(2 * as.matrix(logtheta)[1])          # noise variance
  2 * s2 * diag(dim(x)[1])
}

#' Compute the covariance matrix of a covariance function as the sum of two
#' other covariance functions.
#' @param covFunc1 list with $deriv and $numparams elements impl. a cov. func.
#' @param covFunc2 list with $deriv and $numparams elements impl. a cov. func.
#' @param logtheta vector of log hyperparameters.
#' @param x matrix of training inputs of size n x D.
#' @return covariance matrix.
cov.gram.sum <- function(covFunc1, covFunc2, logtheta,  x) {
  j <- c(covFunc1$numParams(x), covFunc2$numParams(x))  # number of parameters
  v <- NULL
  for (i in 1:length(j)){
    v <- cbind(v, array(rep(i,j[i]), dim=c(1,j[i])))
  }
  A1 <- covFunc1$gram(logtheta[v==1], x)
  A2 <- covFunc2$gram(logtheta[v==2], x)

  A1 + A2
}

#' Compute the derivative matrix of a covariance function as the sum of two
#' other covariance functions.
#' @param covFunc1 list with $deriv and $numparams elements impl. a cov. func.
#' @param covFunc2 list with $deriv and $numparams elements impl. a cov. func.
#' @param logtheta vector of log hyperparameters.
#' @param x matrix of training inputs of size n x D.
#' @param z index of hyperparam in logtheta to compute deriv. with respect to.
#' @return derivative matrix.
cov.deriv.sum <- function(covFunc1, covFunc2, logtheta, x, z) {
  j <- c(covFunc1$numParams(x), covFunc2$numParams(x))  # number of parameters
  v <- NULL
  for (i in 1:length(j)) {
    v <- cbind(v, array(rep(i,j[i]), dim=c(1,j[i])))
  }
  i <- v[z]
  j <- sum(v[1:z] == i)

  if (i == 1) {
    covFunc1$deriv(logtheta[v == i], x, j)
  } else {
    covFunc2$deriv(logtheta[v == i], x, j)
  }
}

#' Computes pairwise squared distances between vectors stored in columns
#' of two matrices.
#' @param A matrix of size D x n.
#' @param B matrix of size D x m.
#' @return matrix of size n x m.
sq_dist <- function(A, B) {
  if(dim(A)[1] != dim(B)[1]){
    stop("Error: column lengths must agree in sq_dist")
  }

  if(length(A) == 1 && length(B) == 1){
    A = as.vector(A)
    B = as.vector(B)
  }

  n <- dim(A)[2]
  m <- dim(B)[2]
  C <- array(0, dim=c(n,m))

  # for m == 1 R automatically turns a column matrix into a row matrix
  if (m == 1) {
    for(d in 1:dim(A)[1]){
      C <- C + t((B[rep(d,n),] - (t(A[rep(d,m),])))^2)
    }
  } else {
    for(d in 1:dim(A)[1]) {
      C <- C + (B[rep(d,n),] - (t(A[rep(d,m),])))^2
    }
  }

  C
}
