\name{simData}
\alias{simData}
\docType{data}
\title{Example dataset for tests
}
\description{
Example dataset for tests
}
\usage{data("simData")}
\format{
  A list with the following entries
  \itemize{
   \item \code{X} Dataframe with 500 observations and three variables.
   \item \code{environment} A vector of length 500, indicating which environment
   the observations belong to.
   \item \code{interventionVar} A matrix of dimension 6 (no. of environments) x 3 (no. of variables),
    where entry i,j indicates whether variable j was intervened on in environment i.
  }
}
\keyword{datasets}
