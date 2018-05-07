library(CondIndTests)
require(testthat)
context("All supported methods")

## TODO: add different subroutines
verbose <- FALSE
dat_gens <- c(1, 2)
methods <- c(
             "KCI",
             "InvariantConditionalQuantilePrediction",
             "InvariantEnvironmentPrediction",
             "InvariantResidualDistributionTest", 
             "InvariantTargetPrediction",
             "ResidualPredictionTest"
             )

settings <- c("1", "2", "3", "4", "5", "6", "7", "8")

for(dat_gen in dat_gens){
  for(setting in settings){
  
    switch(setting, 
             # all arguments one-dimensional
             # second argument continuous
             "1" = {
               dx <- 1
               dy <- 1
               dz <- 1
               n <- 200
               y_type <- "cont"
             },
             # second argument discrete
             "2" = {
               dx <- 1
               dy <- 1
               dz <- 1
               n <- 200
               y_type <- "discrete"
             },
             # conditioning argument two-dimensional
             # second argument continuous
             "3" = {
               dx <- 1
               dy <- 1
               dz <- 2
               n <- 200
               y_type <- "cont"
             },
             # second argument discrete
             "4" = {
               dx <- 1
               dy <- 1
               dz <- 2
               n <- 200
               y_type <- "discrete"
             },
             # second argument two-dimensional
             # second argument continuous
             "5" = {
               dx <- 1
               dy <- 2
               dz <- 2
               n <- 200
               y_type <- "cont"
             },
             # second argument discrete
             "6" = {
               dx <- 1
               dy <- 2
               dz <- 2
               n <- 200
               y_type <- "discrete"
             },
             # first argument two-dimensional
             # second argument continuous
             "7" = {
               dx <- 2
               dy <- 2
               dz <- 2
               n <- 200
               y_type <- "cont"
             },
             # second argument discrete
             "8" = {
               dx <- 2
               dy <- 2
               dz <- 2
               n <- 200
               y_type <- "discrete"
             }
           )
  
    if(dat_gen == 1){
      X <- matrix(rnorm(n*dx), n, dx)
      if(y_type == "cont"){
        Y <- matrix(rnorm(n*dy), n, dy)
      }else{
        Y <- data.frame(matrix(rbinom(n*dy, size = 1, prob = 0.5), n, dy))
        Y <- data.frame(lapply(Y, as.factor))
      }
      Z <- matrix(rnorm(n*dz), n, dz)
    }else{
      X <- if(dx == 1) rnorm(n) else matrix(rnorm(n*dx), n, dx)
      if(y_type == "cont"){
        Y <- if(dy == 1) rnorm(n) else matrix(rnorm(n*dy), n, dy)
      }else{
        if(dy == 1){
          Y <- as.factor(rbinom(n, size = 1, prob = 0.5))
        }else{
          Y <- data.frame(matrix(rbinom(n*dy, size = 1, prob = 0.5), n, dy))
          Y <- data.frame(lapply(Y, as.factor))
        }
      }
      Z <-  if(dz == 1) rnorm(n) else matrix(rnorm(n*dz), n, dz)
    }
    
    
    for(method in methods){
      
      if(method == "KCI" & setting %in% c("1", "2", "3", "4", "5", "6", "7", "8") |
         method == "InvariantConditionalQuantilePrediction" & setting %in% c("2", "4", "6", "8") |
         method == "InvariantEnvironmentPrediction" & setting %in% c("2", "4", "6", "8")  |
         method == "InvariantResidualDistributionTest" & setting %in% c("2", "4", "6", "8") |
         method == "InvariantTargetPrediction" & setting %in% c("1", "2", "3", "4", "5", "6", "7", "8") |
         method == "ResidualPredictionTest" & setting %in% c("1", "2", "3", "4", "5", "6", "7", "8")
         ){
  
        test_that(paste("Checks output type for", method, "and setting", setting), {
          cat(paste("\nMethod:", method, "; Setting", setting, "\n"))
          expect_is(
            out <- CondIndTest(X, Y, Z, method = method)
            , "list")
          if(verbose) print(out)
        }
        )
      }
      
      if(method == "InvariantConditionalQuantilePrediction" & setting %in% c("1", "3", "5", "7") |
         method == "InvariantEnvironmentPrediction" & setting %in% c("1", "3", "5", "7") |
         method == "InvariantResidualDistributionTest" & setting %in% c("1", "3", "5", "7") 
         ){
        
        test_that(paste("Checks output type for", method, "and setting", setting), {
        cat(paste("\nMethod:", method, "; Setting", setting, "\n"))
          expect_error(
        out <- CondIndTest(X, Y, Z, method = method)
        )
        }
        )
      }
      
    }
  }
}