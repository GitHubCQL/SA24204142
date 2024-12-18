% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{GibbsC}
\alias{GibbsC}
\title{A Gibbs sampler using Rcpp}
\usage{
GibbsC(N, a, b, n)
}
\arguments{
\item{N}{the number of samples}

\item{a}{parameter of target density}

\item{b}{parameter of target density}

\item{n}{parameter of target density}
}
\value{
a random sample of size \code{N}
}
\description{
A Gibbs sampler using Rcpp for HW9
}
\examples{
\dontrun{
rnC <- GibbsC(1000, 1, 1, 10)
par(mfrow=c(2,1));
plot(rnC[,1],type='l')
plot(rnC[,2],type='l')
}
}
