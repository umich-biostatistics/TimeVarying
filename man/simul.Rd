% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simul.R
\name{simul}
\alias{simul}
\title{Simulate Data}
\usage{
simul(N = 1000, N_Strata = 10, p = 5)
}
\arguments{
\item{N}{the sample size for each strata.}

\item{N_Strata}{number of stratum.}

\item{p}{number of parameters.}

\item{p_true}{number of true parameters.}
}
\value{
return a list of the following values:
\itemize{
  \item delta:  event indicator
  \item z: Covariate matrix
  \item facility:
  \item time: the death time
}
}
\description{
This simulation function is used to simulate data for testing 
the functions in this package.
}
\examples{
# generate the simuluation data
data2 <- simul(N = 2000, N_Strata = 5, p=5 )

}
