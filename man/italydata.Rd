\name{italydata}
\alias{italydata}
\docType{data}
\title{
COVID-19 Mortality Rates for Italy
}
\description{
A data set containing COVID-19 mortality rates for Italy for a period of 59 days from 27 Feb 2020 to 27 April 2020.
}
\usage{data("italydata")}
\format{
  A data frame with 59 observations on the following 2 variables.
  \describe{
    \item{\code{date}}{a character vector}
    \item{\code{rate}}{a numeric vector}
  }
}

\source{
https://covid19.who.int/
}

\examples{
data(italydata)
## maybe str(italydata) ; plot(italydata) ...
}
\keyword{datasets}
