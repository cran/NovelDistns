\name{ukdata}
\alias{ukdata}
\docType{data}
\title{
COVID-19 Mortality Rates for United Kingdom
}
\description{
A data set containing COVID-19 mortality rates for United Kingdom for a period of 76 days from 15 April 2020 to 30 June 2020
}
\usage{data("ukdata")}
\format{
  A data frame with 76 observations on the following 2 variables.
  \describe{
    \item{\code{date}}{a character vector}
    \item{\code{rate}}{a numeric vector}
  }
}

\source{
https://covid19.who.int/
}

\examples{
data(ukdata)
## maybe str(ukdata) ; plot(ukdata) ...
}
\keyword{datasets}
