\name{kenyadata}
\alias{kenyadata}
\docType{data}
\title{
COVID-19 daily cases for Kenya
}
\description{
A data set containing COVID-19 daily cases for Kenya for a period of 56 days from 28 March 2020 to 24 May 2020
}
\usage{data("kenyadata")}
\format{
  A data frame with 58 observations on the following 2 variables.
  \describe{
    \item{\code{date}}{a character vector}
    \item{\code{cases}}{a numeric vector}
  }
}

\source{
https://covid19.who.int/
}

\examples{
data(kenyadata)
## maybe str(kenyadata) ; plot(kenyadata) ...
}
\keyword{datasets}
