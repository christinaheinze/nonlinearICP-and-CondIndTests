library(CondIndTests)
# context("All supported methods")

##
dat_gen <- 2
methods <- c(
             # "KCI", 
             # "InvariantConditionalQuantilePrediction",
             # "InvariantEnvironmentPrediction",
             # "InvariantResidualDistributionTest",
             # "InvariantTargetPrediction",
             "ResidualPredictionTest")

settings <- c("1", "2", "3", "4", "5", "6", "7", "8")

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
    
    if(setting > 0 & method %in% c("KCI", 
                                    "InvariantConditionalQuantilePrediction",
                                    "InvariantEnvironmentPrediction",
                                    "InvariantResidualDistributionTest",
                                    "InvariantTargetPrediction",
                                    "ResidualPredictionTest")){
      
      # test_that(paste("Checks output type for", method, "and setting", setting), {
        cat(paste("\nMethod:", method, "; Setting", setting, "\n"))
        # expect_is(
          out <- CondIndTest(X, Y, Z, method = method)
          # , "list")
          print(out)
      # }
      # )
    }
  }
}