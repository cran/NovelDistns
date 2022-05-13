\name{bladderdata}
\alias{bladderdata}
\docType{data}
\title{
Bladder Cancer data
}
\description{
A data set containing remission time in months of a sample of 128 bladder cancer patients
}
\usage{data("bladderdata")}
\format{
  A data frame with 128 observations on the following variable.
  \describe{
    \item{\code{time}}{a numeric vector}
  }
}

\source{
E. T. Lee and J. Wang, Statistical Methods for Survival Data Analysis, vol. 476, John Wiley & Sons, Hoboken, NJ, USA, 2003.
}

\examples{
data(bladderdata)
## maybe str(bladderdata) ; plot(bladderdata) ...
}
\keyword{datasets}
