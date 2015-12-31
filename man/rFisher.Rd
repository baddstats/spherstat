\name{rFisher}
\alias{rFisher}
\title{
Generate Fisher distribution on the sphere
}
\description{
Generates variates of the Fisher distribution on the sphere. 
}
\usage{
rFisher(n, mode, kappa, win=sphwin(type="sphere"), nsim=1, drop=TRUE, as.sp=TRUE, ndim="2")
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
  \item{nsim}{
    Number of simulated realisations to be generated.
  }
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the result will be a point pattern, rather than a list containing a point pattern. 
  }
  \item{as.sp}{
    Logical.  If TRUE, returns an object of class as defined by
    \code{sp.dim}.  Otherwise, returns a matrix.  See Value. 
  }
  \item{ndim}{
    A string, taking value \code{"2"} or \code{"3"}.  Specifies whether the object should contain the locations of the points in spherical coordinates (\code{ndim="2"}) or Cartesian coordinates (\code{ndim="3"}).
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
  If \code{nsim=1} and \code{drop=FALSE} then a single item as described below; otherwise a list containing \code{nsim} items.

  An item is determined by the values of \code{as.sp} and \code{ndim}:
  If \code{as.sp=FALSE} and \code{ndim="2"}, a two column matrix giving the locations of the simulated points.
  If \code{as.sp=FALSE} and \code{ndim="3"}, a three column matrix giving the locations of the simulated points.
  If \code{as.sp=TRUE} and \code{ndim="2"}, an object of class \code{\link{sp2}} giving the locations of the simulated points.
  If \code{as.sp=TRUE} and \code{ndim="3"}, an object of class \code{\link{sp3}} giving the locations of the simulated points.
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

rFis2 <- rFisher(n=100, mode=c(pi/3, pi/3), kappa=0.1, win=sphwin(type="band", param=c(0, pi/3), ref=c(pi/3, pi/3)))
rFis2
}
\keyword{spatial}
\keyword{datagen}

