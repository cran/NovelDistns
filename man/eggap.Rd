\name{eggap}
\alias{reggap} \alias{qeggap} \alias{peggap} \alias{deggap} \alias{mleggap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exponentiated Generalized Gull Alpha Power Family of distribution
}
\description{
Computes the pdf, cdf, quantile, and random numbers and estimates the parameters of the exponentiated  \code{G} gull alpha power family of distribution due to Kilai et al. (2022) specified by the cdf.
\deqn{F(x,{\Theta}) = \left[1-\left(1-\frac{\alpha G(x)}{\alpha^{G(x)}}\right)^{a}\right]^{b}} where \eqn{\theta} is the baseline family parameter vector. Also, a>0, b>0 are the extra parameters induced to the baseline cumulative distribution function (cdf) \code{G} whose pdf is \code{g}.
Here, the baseline \code{G} refers to the cdf of: exponential, rayleigh and weibull.
}
\usage{
reggap(n, dist, param)
qeggap(p, dist, param, log.p = FALSE, lower.tail = TRUE)
peggap(data, dist, param, log.p = FALSE, lower.tail = TRUE)
deggap(data, dist, param, log = FALSE)
mleggap(data, dist,starts, method="SANN")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{n}{number of realizations to be generated.}
\item{p}{quantile value between 0 and 1.}
\item{data}{Vector of observations.}
\item{param}{parameter vector \eqn{\Theta=(a,b,\theta,\alpha)}}
\item{log}{If \code{TRUE}, then log(pdf) is returned.}
\item{log.p}{If \code{TRUE}, then log(cdf) is returned and quantile is computed for \code{exp(-p)}.}
\item{lower.tail}{If \code{FALSE}, then \code{1-cdf} is returned and quantile is computed for \code{1-p}.}
\item{dist}{The name of family's pdf including: "\code{exponential}", "\code{rayleigh}", "\code{weibull}", "\code{lomax}"}
\item{method}{the method for optimizing the log likelihood function.  It can be one of \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"} or \code{"SANN"}.  The default is \code{"BFGS"}.  The details of these methods can be found in the manual pages for  \code{optim}}
\item{starts}{initial values of \code{(theta, a, b, alpha)}}
}

\value{\enumerate{
\item A vector of the same length as \code{data}, giving the pdf values computed at \code{data}.
\item A vector of the same length as \code{data}, giving the cdf values computed at \code{data}.
\item A vector of the same length as \code{p}, giving the quantile values computed at \code{p}.
\item A vector of the same length as \code{n}, giving the random numbers realizations.
\item A sequence of goodness-of-fit statistics such as: Akaike Information Criterion (\code{AIC}), Consistent Akaike Information Criterion (\code{CAIC}), Bayesian Information Criterion (\code{BIC}), Hannan-Quinn information criterion (\code{HQIC}), Cramer-von Misses statistic (\code{CM}), Anderson Darling statistic (\code{AD}), log-likelihood statistic (\code{log}). The Kolmogorov-Smirnov (\code{KS}) test statistic and corresponding \code{p-value} and the convergence status.
}}
\references{
Mutua Kilai et al (2022) A new generalization of Gull Alpha Power Family of distributions with application to modeling COVID-19 mortality rates, https://doi.org/10.1016/j.rinp.2022.105339.
}
\author{
Mutua Kilai, Gichuhi A. Waititu, Wanjoya A. Kibira
}


\examples{
x=runif(10,min=0,max=1)
reggap(10,"exp",c(0.3,0.5,0.7,0.8))
qeggap(0.6,"exp",c(0.3,0.5,0.7,0.8))
peggap(x,"exp",c(0.3,0.5,0.7,0.8))
deggap(x,"exp",c(0.3,0.5,0.7,0.8))
mleggap(x,"exp",c(0.3,0.5,0.7,0.8))
}
