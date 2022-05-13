\name{egap}
\alias{regap} \alias{qegap}\alias{pegap} \alias{degap} \alias{mlegap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exponentiated Gull Alpha Power Family of distribution
}
\description{
Computes the pdf, cdf, quantile, and random numbers and estimates the parameters of the exponentiated  gull alpha power family of distribution specified by the cdf.
\deqn{F(x,{\Theta}) = \left[\frac{\alpha G(x)}{\alpha^{G(x)}}\right]^{b}} where \eqn{\theta} is the baseline family parameter vector. Also, b>0 are the extra parameters induced to the baseline cumulative distribution function (cdf) \code{G} whose pdf is \code{g}.
Here, the baseline \code{G} refers to the cdf of: exponential, rayleigh and weibull.
}
\usage{
regap(n, dist, param)
qegap(p, dist, param, log.p = FALSE, lower.tail = TRUE)
pegap(data, dist, param, log.p = FALSE, lower.tail = TRUE)
degap(data, dist, param, log = FALSE)
mlegap(data, dist,starts, method="SANN")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{n}{number of realizations to be generated.}
\item{p}{quantile value between 0 and 1.}
\item{data}{Vector of observations.}
\item{param}{parameter vector \eqn{\Theta=(b,\theta,\alpha)}}
\item{log}{If \code{TRUE}, then log(pdf) is returned.}
\item{log.p}{If \code{TRUE}, then log(cdf) is returned and quantile is computed for \code{exp(-p)}.}
\item{lower.tail}{If \code{FALSE}, then \code{1-cdf} is returned and quantile is computed for \code{1-p}.}
\item{dist}{The name of family's pdf including: "\code{exponential}", "\code{rayleigh}", "\code{weibull}", "\code{lomax}"}
\item{method}{the method for optimizing the log likelihood function.  It can be one of \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"} or \code{"SANN"}.  The default is \code{"BFGS"}.  The details of these methods can be found in the manual pages for  \code{optim}}
\item{starts}{initial values of \code{(theta, b, alpha)}}

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

\author{
Mutua Kilai, Gichuhi A. Waititu, Wanjoya A. Kibira
}


\examples{
x=runif(10,min=0,max=1)
regap(10,"exp",c(0.3,0.5,0.7))
qegap(0.6,"exp",c(0.3,0.5,0.7))
pegap(x,"exp",c(0.3,0.5,0.7))
degap(x,"exp",c(0.3,0.5,0.7))
mlegap(x,"exp",c(0.3,0.5,0.7))
}
