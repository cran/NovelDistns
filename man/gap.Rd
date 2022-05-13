\name{gap}
\alias{rgap} \alias{qgap} \alias{pgap} \alias{dgap} \alias{mlgap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gull Alpha Power Family of distribution
}
\description{
Computes the pdf, cdf, quantile, and random numbers and estimates the parameters of the exponentiated  gull alpha power family of distribution specified by the cdf.
\deqn{F(x,{\Theta}) = \left[\frac{\alpha G(x)}{\alpha^{G(x)}}\right]} where \eqn{\theta} is the baseline family parameter vector.Here, the baseline \code{G} refers to the cdf of: exponential, rayleigh and weibull.
}
\usage{
rgap(n, dist, param)
qgap(p, dist, param, log.p = FALSE, lower.tail = TRUE)
pgap(data, dist, param, log.p = FALSE, lower.tail = TRUE)
dgap(data, dist, param, log = FALSE)
mlgap(data, dist,starts, method="SANN")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{n}{number of realizations to be generated.}
\item{p}{quantile value between 0 and 1.}
\item{data}{Vector of observations.}
\item{param}{parameter vector \eqn{\Theta=(\theta,\alpha)}}
\item{log}{If \code{TRUE}, then log(pdf) is returned.}
\item{log.p}{If \code{TRUE}, then log(cdf) is returned and quantile is computed for \code{exp(-p)}.}
\item{lower.tail}{If \code{FALSE}, then \code{1-cdf} is returned and quantile is computed for \code{1-p}.}
\item{dist}{The name of family's pdf including: "\code{exponential}", "\code{rayleigh}", "\code{weibull}", "\code{lomax}"}
\item{method}{the method for optimizing the log likelihood function.  It can be one of \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"} or \code{"SANN"}.  The default is \code{"BFGS"}.  The details of these methods can be found in the manual pages for  \code{optim}}
\item{starts}{initial values of \code{(theta, alpha)}}
}

\value{
\enumerate{
\item A vector of the same length as \code{data}, giving the pdf values computed at \code{data}.
\item A vector of the same length as \code{data}, giving the cdf values computed at \code{data}.
\item A vector of the same length as \code{p}, giving the quantile values computed at \code{p}.
\item A vector of the same length as \code{n}, giving the random numbers realizations.
\item A sequence of goodness-of-fit statistics such as: Akaike Information Criterion (\code{AIC}), Consistent Akaike Information Criterion (\code{CAIC}), Bayesian Information Criterion (\code{BIC}), Hannan-Quinn information criterion (\code{HQIC}), Cramer-von Misses statistic (\code{CM}), Anderson Darling statistic (\code{AD}), log-likelihood statistic (\code{log}). The Kolmogorov-Smirnov (\code{KS}) test statistic and corresponding \code{p-value} and the convergence status.
}
}
\references{
Muhammad et al (2020) A Gull Alpha Power Weibull distribution with applications to real and simulated data. https://doi.org/10.1371/journal.pone.0233080
}
\author{
Mutua Kilai, Gichuhi A. Waititu, Wanjoya A. Kibira
}

\examples{
x=runif(10,min=0,max=1)
rgap(10,"exp",c(0.3,0.5))
qgap(0.6,"exp",c(0.3,0.5))
pgap(x,"exp",c(0.3,0.5))
dgap(x,"exp",c(0.3,0.5))
mlgap(x,"exp",c(0.3,0.5))
}
