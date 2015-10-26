\name{rFisher}
\alias{rFisher}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Fisher distribution on the sphere
}
\description{
Generates variates of the Fisher distribution on the sphere. 
}
\usage{
rFisher(n, mode, kappa, win=sphwin())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
The number of points to be simulated.  An integer.
}
  \item{mode}{
A length 2 numeric giving the mode of the distribution
}
  \item{kappa}{
The concentration parameter.  Must be positive.
}
  \item{win}{
The window, an object of type sphwin.
}
}

\details{
\code{rFisher} generates a realisation of \code{n} points of the Fisher distribution inside the window \code{win}.  The Fisher distribution is the analogue on the sphere of the normal distribution; the maximum of the density occurs at \code{mode} and the minimum at the point diametrically opposite \code{mode} (i.e. the antimode).  The density of this distribution on any small or great circle that is orthogonal to the axis containing the mode and antimode, is uniform.
}
\value{
An object of class \code{\link{sp2}}, giving the locations of the simulated points.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
rFis <- rFisher(n=100, mode=c(0,0), kappa=50, win=sphwin())
rFis

rFis2 <- rFisher(n=100, mode=c(pi/3, pi/3), kappa=0.1, win=sphwin(type="cap", param=c(0, pi/3), ref=c(pi/3, pi/3)))
rFis2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
