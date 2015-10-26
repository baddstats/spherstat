\name{rdiscunif}
\alias{rdiscunif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The discrete uniform distribution
}
\description{
Generate variates from a discrete uniform distribution
}
\usage{
rdiscunif(n, min, max, steps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
The number of points to be simulated.  An integer.
}
  \item{min}{
The minimum value of the support of the distribution.
}
  \item{max}{
The maximum value of the support of the distribution.
}
  \item{steps}{
The step size in the support of the distribution i.e. the difference between consecutive values in the support.  Must be positive and a factor of \code{(max-min)}.
}
}

\details{
\code{rdiscunif} generates a realisation of \code{n} points of the discrete uniform distribution with support \code{(min, min+steps, min+2*steps, ... max)}.  If \code{steps} is not a factor of \code{(max-min)}; or if \code{min > max} then an error is generated.  If \code{min = max}, then \code{rdiscunif(n, min, max, steps)=repeat(min, n)}
}

\value{
A length \code{n} numeric containing realisations of the specified uniform distribution.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tom Lawrence tjlawrence@bigpond.com
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{runif}} for generating the continuous uniform distribution.
}
\examples{
rdiscunif(n=10, min=0, max=5, steps=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
