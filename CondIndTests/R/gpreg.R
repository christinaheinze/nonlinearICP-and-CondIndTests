#' R implemantation of hyperparameter optimization in Gaussian Process (GP)
#' regression based on vers. 2 of the GPML toolbox by Carl Edward Rasmussen and
#' Christopher K. I. Williams distributed with their book "Gaussian Processes
#' for Machine Learning", the MIT Press 2006.
#'
#' Copyright (c) 2005-2018 Carl Edward Rasmussen & Hannes Nickisch. All rights
#' reserved.
#'
#' Original code was distibuted under the FreeBSD license. See
#' http://www.gaussianprocess.org/gpml/code/matlab/Copyright for additional
#' information.
#'

#' Find optimal hyperparameters by maximizing the marginal likelihood.
#' @param start_logtheta initial log hyperperameter values.
#' @param x matrix of training inputs of size n x D.
#' @param y vector of targets of size n.
#' @param covfunc covariance function defining the gaussian process.
#' @param verbose If TRUE, log information about the progress of the
#'                optimization to the console
#' @return optimal hyperparameter values.
gpr.find_optimal_hyperparams <- function(
  start_logtheta, x, y, covfunc, verbose = FALSE
) {
  res <- mize(
    start_logtheta,                 # initial hyperparameter values
    gpr.optim_funcs(x, y, covfunc), # closure of functions for optimization
    method = 'CG',                  # use Conjugate Gradients method
    verbose = verbose               # print optimization output?
  )

  res$par                           # return optimal parameter values
}

#' Closure of covariance function and derivative of covariance function to be
#' used by optimization function.
#' @param x matrix of training inputs of size n x D.
#' @param y vector of targets of size n.
#' @param covfunc covariance function defining the gaussian process.
#' @return list of closures.
gpr.optim_funcs <- function(x, y, covfunc) {
  list(
    # compute negative log marginal likelihood (function to be optimized)
    fn = function(logtheta) {
      C <- inv_cov_matrix(logtheta, x, y, covfunc)
      neg_log_marginal_likelihood(C$alpha, y, C$L)
    },

    # compute gradient of negative log marginal likelihood
    gr = function(logtheta) {
      C <- inv_cov_matrix(logtheta, x, y, covfunc)
      grad_neg_log_marginal_likelihood(logtheta, C$alpha, C$L, x, covfunc)
    },

    # compute both function and gradient of negative log marginal likelihood
    # (reusing the same matrix inversion for both)
    fg = function(logtheta) {
      C <- inv_cov_matrix(logtheta, x, y, covfunc)
      list(
        fn = neg_log_marginal_likelihood(C$alpha, y, C$L),
        gr = grad_neg_log_marginal_likelihood(logtheta, C$alpha, C$L, x, covfunc)
      )
    }
  )
}

#' Invert covariance matrix using cholesky decomposition.
#' @param logtheta log hyperperameter values.
#' @param x matrix of training inputs of size n x D.
#' @param y vector of targets of size n.
#' @param covfunc covariance function defining the gaussian process.
inv_cov_matrix <- function(logtheta, x, y, covfunc) {
  logtheta <- as.matrix(logtheta)

  K <- covfunc$gram(logtheta, x)       # training set covariance matrix
  L <- t(cholesky(K))                  # cholesky factorization of covariance
  alpha <- solve_cholesky(t(L), y)

  list(L = L, alpha = alpha)
}

#' Compute negative log marginal likelihood.
#' @param alpha solution of the linear system of equations \eqn{\alpha = L^T (L y)}.
#' @param y vector of targets of size n.
#' @param L Cholesky factor of training set covariance matrix.
neg_log_marginal_likelihood <- function(alpha, y, L) {
  n <- dim(y)[1]
  m <- dim(y)[2]

  0.5 * Trace(t(y) %*% alpha) + m * sum(log(diag(L))) + 0.5 * m * n * log(2*pi)
}

#' Compute gradient of negative log marginal likelihood.
#' @param logtheta log hyperperameter values.
#' @param alpha solution of the linear system of equations \eqn{\alpha = L^T (L y)}.
#' @param L Cholesky factor of training set covariance matrix.
#' @param x matrix of training inputs of size n x D.
#' @param covfunc covariance function defining the gaussian process.
grad_neg_log_marginal_likelihood <- function(logtheta, alpha, L, x, covfunc) {
  dnlml <- matrix(0, dim(as.matrix(logtheta))[1])

  n <- dim(x)[1]
  m <- dim(alpha)[2]
  W <- m * solve(t(L), solve(L, eye(n))) - alpha %*% t(alpha)

  for (i in 1:length(dnlml)) {
    dK <- covfunc$deriv(logtheta, x, z = i)
    dnlml[i] <-  sum(W * dK)/2
  }

  as.vector(dnlml)
}

#' Compute the Choleski factorization of a real symmetric positive-definite
#' square matrix using R's `chol`.
#' @param K numeric (or logical) symmetric, positive-definite matrix.
#' @return cholesky factorization of covariance.
cholesky <- function(K) { as.matrix(chol(K)) }

#' Solve sytem of linear equations from the Cholesky factorization,
#' \eqn{\alpha = A^T (A B)}, where A is square, symmetric, positive definite.
#' @param A square, symmetric, positive definite matrix.
#' @param B vector or matrix.
solve_cholesky <- function(A, B) {
  if (dim(A)[1] != dim(A)[2]) {
    stop(sprintf('A needs to be symmetric: (%d, %d)', dim(A)))
  } else if (dim(A)[1] != dim(B)[1]) {
    stop(sprintf(
      'Wrong sizes of matrix arguments A and B; A: (%d, %d); B: (%d, %d);'
    , append(dim(A), dim(B))))
  }

  # solve
  X <- solve(A, solve(t(A), B))

  as.matrix(X)
}
