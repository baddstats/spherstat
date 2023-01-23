\name{rMaternII.sphwin}
\alias{rMaternII.sphwin}
\title{
  Simulate Matern Model II
}
\description{
  Generate a random point pattern, a simulated realisation of the Matern
  Model II inhibition process model.  
}
\usage{
rMaternII.sphwin(kappa, r, win = sphwin(type="sphere"), nsim=1, drop=TRUE,
 stationary=TRUE, as.sp=TRUE, ndim="2")
}
\arguments{
  \item{kappa}{
    Intensity of the Poisson process of proposal points (a positive number).
  }
  \item{r}{
    Inhibition distance
  }
   \item{win}{
     The window in which to simulate the pattern, an object of type
     \code{\link{sphwin}}. 
  }
  \item{nsim}{
    Number of simulated realisations to be generated.
  }
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the result will be a point pattern, rather than a list containing a point pattern. 
  }
  \item{stationary}{ 
    Logical. If \code{TRUE} (the default), the process is simulated on the entire sphere, and the output is only those points in the window defined by \code{win}.  If \code{FALSE}, the process is simulated within the window defined by \code{win}.
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
  This algorithm generates a realisation of Matern's Model II inhibition
  process inside the window \code{win}.  First, a simulated stationary
  Poisson process with intensity \code{kappa} is generated.  The number of
  proposal points will equal the number of points in this pattern; the
  original pattern is ignored after this step. 

  The simulated pattern is comprised of proposal points.  The first
  proposal point is automatically provisionally retained. The remaining
  proposal points are tested one at a time.  If a proposal point occurs
  within great circle distance \code{r} of any existing provisionally
  retained points, those points are deleted and the proposal point
  rejected.  Otherwise, the proposal point is provisionally retained. 

  Once all the proposal points have been tested, the final pattern of
  retained points constitute Matern's Model II. 

  The process is constructed by first generating a uniform Poisson point process of "proposal" points with intensity \code{kappa}. If \code{stationary = TRUE} (the default), the proposal points are generated on the entire sphere; this is relevant when there may be points outside the window, as these would affect whether points inside the window are retained in the final model. If \code{stationary=FALSE} then the proposal points are only generated inside the window \code{win}; this is relevant when it is known or assumed that there are no points outside the window. 

}
\value{
  If \code{nsim=1} and \code{drop=FALSE} then a single item as described below; otherwise a list containing \code{nsim} items.

  An item is determined by the values of \code{as.sp} and \code{ndim}:

  If \code{as.sp=FALSE} and \code{ndim="2"}, a two column matrix giving the locations of the simulated points.

  If \code{as.sp=FALSE} and \code{ndim="3"}, a three column matrix giving the locations of the simulated points.

  If \code{as.sp=TRUE} and \code{ndim="2"}, an object of class \code{sp2} giving the locations of the simulated points.

  If \code{as.sp=TRUE} and \code{ndim="3"}, an object of class \code{sp3} giving the locations of the simulated points.
}
\references{
  Lawrence, T.J. (20178) Master's Thesis, University of Western Australia.
}
\author{
Tom Lawrence <email:tjlawrence@bigpond.com>
}
\note{
This function is the analogue for point processes on the sphere of the
  function \code{\link[spatstat.random]{rMaternII}} in \pkg{spatstat} for point processes in \eqn{R^2}{R^2}.  Hence
  elements of this help page have been taken from that for 
  \code{\link[spatstat.random]{rMaternII}}, with the permission of
  A. J. Baddeley.  This enables the information on this help page to be
  consistent with that for\code{\link[spatstat.random]{rMaternII}}.  It is
  hoped that this will minimise or 
  remove any confusion for users of both \pkg{spatstat} and
  \pkg{spherstat}. 
}
\seealso{
  \code{\link[spatstat.random]{rHardcore}}, \code{\link{rMatClust.sphwin}},
  \code{\link{rMaternI.sphwin}}, \code{\link[spatstat.random]{rMaternII}},
  \code{\link{rpoispp.sphwin}}, \code{\link{rStrauss.sphwin}},
  \code{\link{rThomas.sphwin}} 
}
\examples{
rM2 <- rMaternII.sphwin(kappa=5, r=pi/20, win=sphwin())
rM2
}
\keyword{spatial}
\keyword{datagen}

