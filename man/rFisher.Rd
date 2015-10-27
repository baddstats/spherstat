\name{rFisher}
\alias{rFisher}
\title{
Generate Fisher distribution on the sphere
}
\description{
Generates variates of the Fisher distribution on the sphere. 
}
\usage{
rFisher(n, mode, kappa, win=sphwin())
}
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
  \code{rFisher} generates a realisation of \code{n} points of the
  Fisher distribution inside the window \code{win}.  The Fisher
  distribution is the analogue on the sphere of the normal distribution;
  the maximum of the density occurs at \code{mode} and the minimum at
  the point diametrically opposite \code{mode} (i.e. the antimode).  The
  density of this distribution on any small or great circle that is
  orthogonal to the axis containing the mode and antimode, is uniform. 
}
\value{
  An object of class \code{\link{sp2}}, giving the locations of the
  simulated points. 
}
\references{
  Lawrence, T.J. (2016) Master's Thesis, University of Western Australia.  
}
\author{
Tom Lawrence <email:tjlawrence@bigpond.com>
}
\seealso{
  \code{\link{runif.sphwin}}
}
\examples{
rFis <- rFisher(n=100, mode=c(0,0), kappa=50, win=sphwin())
rFis

rFis2 <- rFisher(n=100, mode=c(pi/3, pi/3), kappa=0.1, win=sphwin(type="cap", param=c(0, pi/3), ref=c(pi/3, pi/3)))
rFis2
}
\keyword{spatial}
\keyword{datagen}

