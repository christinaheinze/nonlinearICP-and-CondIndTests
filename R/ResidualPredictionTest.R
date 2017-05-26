#' Residual prediction test.
#'
#' @description Tests the null hypothesis that Y and E are independent given X.
#'
#' @param Y An n-dimensional vector.
#' @param E An n-dimensional vector or an nxq dimensional matrix or dataframe.
#' @param X A matrix or dataframe with n rows and p columns.
#' @param alpha Significance level. Defaults to 0.05.
#' @param verbose If \code{TRUE}, intermediate output is provided. Defaults to \code{FALSE}.
#' @param test
#' @param degree
#' @param basis
#' @param resid_type
#' @param XBasis
#' @param noiseMat
#' @param getnoiseFct
#' @param argsGetNoiseFct
#' @param nSim
#' @param funcOfRes
#' @param useX
#' @param returnXBasis
#' @param nSub
#' @param ntree Random forest parameter: Number of trees to grow. Defaults to 500.
#' @param nodesize Random forest parameter: Minimum size of terminal nodes. Defaults to 5.
#' @param maxnodes Random forest parameter: Maximum number of terminal nodes trees in the forest can have.
#' Defaults to NULL.
#' @param nSeqTests Bonferroni adjustment factor if previous tests where performed
#' (e.g. with subsamples).
#'
#' @return A list with the following entries:
#' \itemize{
#'  \item \code{pValue} The p value for the null hypothesis that Y and E are independent given X.
#'  \item \code{XBasis}
#'  \item \code{fctBasisExpansion}
#'  }
#'
#' @examples
#' n <- 500
#' E <- rbinom(n, size = 1, prob = 0.2)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * (X)^2 + rnorm(n)
#' ResidualPredictionTest(Y, as.factor(E), X)
#'
#' E <- rbinom(n, size = 1, prob = 0.2)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * E + rnorm(n)
#' ResidualPredictionTest(Y, as.factor(E), X)
#'
#' E <- rnorm(n)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * (X)^2 + rnorm(n)
#' ResidualPredictionTest(Y, E, X)
#' ResidualPredictionTest(Y, X, E)

ResidualPredictionTest <- function(Y, E, X,
                            alpha = 0.05,
                            verbose = FALSE,
                            test = NULL,
                            degree = 4,
                            basis = c("nystrom",
                                      "nystrom_poly",
                                      "fourier",
                                      "polynomial",
                                      "provided")[1],
                            resid_type = "OLS",
                            XBasis = NULL,
                            noiseMat = NULL,
                            getnoiseFct = getNoiseStdNorm,
                            argsGetNoiseFct = NULL,
                            nSim = 100,
                            funcOfRes = function(x){abs(x)},
                            useX = TRUE,
                            returnXBasis = FALSE,
                            nSub = ceiling(NROW(X)/4),
                            ntree = 100,
                            nodesize = 5,
                            maxnodes = NULL,
                            nSeqTests = 1,
                            returnModel = FALSE){

  if(verbose){
    cat("\nFunction of residuals:")
    print(body(funcOfRes))
  }

  n <- NROW(X)
  p <- NCOL(X)
  dimE <- NCOL(E)

  if(verbose){
    cat(paste("\nThe environment variable is", dimE, "-dimensional."))
  }

  if(is.null(noiseMat)){
    noiseMat <- matrix(nrow=nSim,ncol=n)
    for (sim in 1:nSim){
      noiseMat[sim,] <- do.call(getnoiseFct, c(list(n, unlist(argsGetNoiseFct), verbose)))
    }
  }else{
    if(nrow(noiseMat) != nSim){
      stop("The provided noise matrix needs to have nSim rows.")
    }
    if(verbose){
      cat(paste("\nNoise matrix has been provided (",
                nrow(noiseMat), "simulations)."))

    }
  }

  if(NCOL(X) == 1 & all(X == 1)){
    ll <- lm(Y ~ 1)
    res <- residuals(ll)
    res <- res/sd(res)
  }else{
    if(basis != "provided"){
      basisComp <- computeBasis(basis, X, p, degree, nSub, verbose)
      XBasis <- basisComp$XBasis
      fctBasisExpansion <- basisComp$fctBasisExpansion
    }else if(!is.null(XBasis)){
      if(verbose){
        cat(paste("\nBasis has been provided. Design matrix has",
                  ncol(XBasis), "columns."))
      }
    }else{
      stop("\nBasis needs to be provided or estimation procedure needs to be specified.")
    }

    RPtestRes <- try(RPtest(x = XBasis,
                        y = Y,
                        resid_type = resid_type,
                        B = nSim,
                        noise_matrix = t(noiseMat),
                        resid_only = TRUE,
                        verbose = verbose))

    if(inherits(RPtestRes, "try-error")){
      RPtestRes <- try(RPtest(x = XBasis,
                              y = Y,
                              resid_type = "Lasso",
                              B = nSim,
                              noise_matrix = t(noiseMat),
                              resid_only = TRUE,
                              verbose = verbose))

      if(!inherits(RPtestRes, "try-error")){
        cat("\nUsing Lasso residuals as previous computation of residuals failed.")
      }
    }

    res <- RPtestRes$resid
  }

  rf <- try(randomForest(x = if(useX) data.frame(E, X) else data.frame(E),
                         y = res,
                         ntree = ntree,
                         nodesize = nodesize,
                         maxnodes = maxnodes))
  pred <- predict(rf, x = if(useX) data.frame(E, X) else data.frame(E))
  # pred <- rf$predicted

  gRes <- funcOfRes(res)
  rf2 <- try(randomForest(x = if(useX) data.frame(E, X) else data.frame(E),
                          y = gRes,
                          ntree = ntree,
                          nodesize = nodesize,
                          maxnodes = maxnodes))
  pred2 <- predict(rf2, x = if(useX) data.frame(E, X) else data.frame(E))

  # under the null stat should be small: we expect residuals to behave
  # roughly like the noise term, no signal left,
  # so we should not be able to achieve a low prediction error;
  # so 1- prediction error should be small
  # under the alternative we expect it stat
  # to be large: can achieve low prediction error when signal is left in
  # residuals
  stat <- 1-mean( (pred-res)^2)/mean( (res-mean(res))^2)
  stat2 <- 1-mean( (pred2-gRes)^2)/mean( (gRes-mean(gRes))^2)

  statsim <- numeric(nSim)
  statsim2 <- numeric(nSim)


  # for (sim in 1:nSim){
  #
  #   if(NCOL(X) == 1 & all(X == 1)){
  #     simnoise <- noiseMat[sim,]
  #     lsim <- lm(simnoise ~ 1)
  #     res <- residuals(lsim)
  #     res <- res/sd(res)
  #   }else{
  #     res <- RPtestRes$resid_sim[,sim]
  #   }
  #
  #   rf <- randomForest(x = if(useX) data.frame(E, X) else data.frame(E),
  #                      y = res,
  #                      ntree = ntree,
  #                      maxdepth = maxdepth)
  #   pred <- predict(rf, x = if(useX) data.frame(E, X) else data.frame(E))
  #   # pred <- rf$predicted
  #   statsim[sim] <- 1-mean( (pred-res)^2)/mean( (res-mean(res))^2)
  #
  #   gRes <- funcOfRes(res)
  #   rf2 <- randomForest(x = if(useX) data.frame(E, X) else data.frame(E),
  #                       y = gRes,
  #                       ntree = ntree,
  #                       maxdepth = maxdepth)
  #   pred2 <- predict(rf2, x = if(useX) data.frame(E, X) else data.frame(E))
  #   # pred <- rf$predicted
  #   statsim2[sim] <- 1-mean( (pred2-gRes)^2)/mean( (gRes-mean(gRes))^2)
  # }

  inputMat <- if(NCOL(X) == 1 & all(X == 1)) noiseMat else RPtestRes$resid_sim
  ret <- apply(inputMat, 2, function(res){
        if(NCOL(X) == 1 & all(X == 1)){
          simnoise <- res
          lsim <- lm(simnoise ~ 1)
          res <- residuals(lsim)
          res <- res/sd(res)
        }

        rf <- randomForest(x = if(useX) data.frame(E, X) else data.frame(E),
                           y = res,
                           ntree = ntree,
                           maxdepth = maxdepth)
        pred <- predict(rf, x = if(useX) data.frame(E, X) else data.frame(E))
        # pred <- rf$predicted
        statsim1 <- 1-mean( (pred-res)^2)/mean( (res-mean(res))^2)

        gRes <- funcOfRes(res)
        rf2 <- randomForest(x = if(useX) data.frame(E, X) else data.frame(E),
                            y = gRes,
                            ntree = ntree,
                            maxdepth = maxdepth)
        pred2 <- predict(rf2, x = if(useX) data.frame(E, X) else data.frame(E))
        # pred <- rf$predicted
        statsim2 <- 1-mean( (pred2-gRes)^2)/mean( (gRes-mean(gRes))^2)
        c(statsim1, statsim2)
  })

  statsim <- ret[1,]
  statsim2 <- ret[2,]

  pvalue1 <- min(1,mean(statsim>stat) +1/length(statsim))
  pvalue2 <- min(1,mean(statsim2>stat2) +1/length(statsim2))
  pvalue <- min(2*min(pvalue1, pvalue2)*nSeqTests, 1)

  if(verbose){
    cat(paste("\np-value: ", signif(pvalue,2)))
  }

  list(pvalue = pvalue, XBasis = if(returnXBasis) XBasis else NULL,
       fctBasisExpansion =
         if(exists("fctBasisExpansion")) fctBasisExpansion else NULL)
}

getNoiseStdNorm <- function(n,args,verbose = FALSE){
  x <- rnorm(n)
  return(x)
}

