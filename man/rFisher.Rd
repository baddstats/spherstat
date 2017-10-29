\name{rFisher}
\alias{rFisher}
\alias{rWatson}
\alias{rWatson.bipolar}
\alias{rWatson.girdle}
\title{
Generate Fisher and Watson distributions on the sphere
}
\description{
Generates variates of the Fisher and Watson distributions on the sphere. 
}
\usage{
rFisher(n, mode, kappa, win=sphwin(type="sphere"), nsim=1, drop=TRUE, as.sp=TRUE, ndim="2")
rWatson(n, mode, kappa, win=sphwin(type="sphere"), squeeze=FALSE, nsim=1, drop=TRUE, as.sp=TRUE, ndim="2")
rWatson.bipolar(n, mode, kappa, win=sphwin(type="sphere"), squeeze=FALSE, nsim=1, drop=TRUE, as.sp=TRUE, ndim="2")
rWatson.girdle(n, mode, kappa, win=sphwin(type="sphere"), squeeze=FALSE, nsim=1, drop=TRUE, as.sp=TRUE, ndim="2")
}
\arguments{
  \item{n}{
    The number of points to be simulated.  An integer.
  }
  \item{mode}{
    A length 2 numeric giving the mode of the distribution
  }
  \item{kappa}{
    The concentration parameter.  Must be positive for the Fisher and Watson bipolar distributions.  Must be negative for the Watson girdle distribution.
  }
  \item{win}{
    The window, an object of type sphwin.
  }
  \item{squeeze}{
    Logical.  For Watson distributions; if TRUE then applies a squeeze test which can speed up the simulation if \code{kappa} is large.
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
  orthogonal to the axis containing the mode and antimode, is uniform.  The algorithm was given in Fisher et al. (1981) and in Fisher et al. (1993).

  \code{rWatson.bipolar} generates a realisation of of \code{n} points of
  the Watson bipolar distribution inside the window \code{n}.  This 
  distribution is the bipolar analogue of the Fisher distribution, with
  maximum density at two diametrically opposite modes (one of which is
  specified as \code{mode} in the function) and minimum density at the great
  circle halfway between the modes.  The algorithm was given in Best and Fisher (1986) and in Fisher et al. (1993).

  \code{rWatson.girdle} generates a realisation of of \code{n} points of
  the Watson bipolar distribution inside the window \code{n}.  This 
  distribution is the `opposite' of a the Watson bipolar in that the
  minimum density is at two diametrically opposite antimodes (one of which is
  specified as \code{mode} in the function) and maximum density at the great
  circle halfway between the antimodes.  The algorithm was given in Best and Fisher (1986) and in Fisher et al. (1993).

  \code{rWatson} simulates either a Watson bipolar distribution (using \code{rWatson.bipolar})
  or a Watson girdle distribution (using \code{rWatson.girdle}) depending on what value of \code{kappa} is provided.

  The algorithsm used in these simulations are found in Section 3.6.2. of N. I. Fisher et al.
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
  Best, D. J. and Fisher, N. I. (1986) `Goodness-of-fit and discordancy tests for samples from the Watson distribution on the sphere', \emph{Australasian Journal of Statistics}, \bold{28}, 13--31.

  Fisher, N. I., Embleton, B. J. J. and Lewis, T. (1987) `Statistical analysis of spherical data', 1st paperback edition.  Cambridge University Press.

  Fisher, N. I., Lewis, T. and Willcox, M. E. (1981) `Tests of discordancy for samples from Fisher's distribution on the sphere' \emph{Journal of the Royal Statistical Society Series C}, \bold{30}, 230--237.

  Lawrence, T.J. (2017) Master's Thesis, University of Western Australia.  
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

rWB <- rWatson.bipolar(n=100, mode=c(0,0), kappa=50, win=sphwin())
rWB

rWB2 <- rWatson.bipolar(n=100, mode=c(pi/3, pi/3), kappa=0.1, win=sphwin(type="band", param=c(0, pi/3), ref=c(pi/3, pi/3)))
rWB2

rWG <- rWatson.girdle(n=100, mode=c(0,0), kappa=-50, win=sphwin())
rWG

rWG2 <- rWatson.girdle(n=100, mode=c(pi/3, pi/3), kappa=-0.1, win=sphwin(type="band", param=c(0, pi/3), ref=c(pi/3, pi/3)))
rWG2

rWat <- rWatson(n=100, mode=c(0,0), kappa=50, win=sphwin())
rWat

rWat2 <- rWatson(n=100, mode=c(pi/3, pi/3), kappa=-0.1, win=sphwin(type="band", param=c(0, pi/3), ref=c(pi/3, pi/3)))
rWat2

rWatFish <- rWatson(n=100, mode=c(0,0), kappa=0, win=sphwin())
rWatFish

}
\keyword{spatial}
\keyword{datagen}

