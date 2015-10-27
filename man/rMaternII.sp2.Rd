\name{rMaternII.sp2}
\alias{rMaternII.sp2}
\title{
  Simulate Matern Model II
}
\description{
  Generate a random point pattern, a simulated realisation of the Matern
  Model II inhibition process model.  
}
\usage{
rMaternII.sp2(kappa, r, win = sphwin(type="sphere"))
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
}
\value{
  The simulated point pattern (an object of class \code{\link{sp2}}).
}
\references{
  Lawrence, T.J. (2016) Master's Thesis, University of Western Australia.
}
\author{
Tom Lawrence <email:tjlawrence@bigpond.com>
}
\note{
This function is the analogue for point processes on the sphere of the
  function \code{\link[spatstat]{rMaternII}} in \pkg{spatstat}.  Hence
  elements of this help page have been taken from that for
  \code{\link[spatstat]{rMaternII}}, with the permission of
  A. J. Baddeley.  This enables the information on this help page to be
  consistent with that for\code{\link[spatstat]{rMaternII}}.  It is
  hoped that this will minimise or 
  remove any confusion for users of both \pkg{spatstat} and
  \pkg{spherstat}. 
}
\seealso{
  \code{link[spatstat]{rHardcore}}, \code{\link{rMatClust.sp2}},
  \code{\link{rMaternI.sp2}}, \code{link[spatstat]{rMaternII}},
  \code{\link{rpoispp.sp2}}, \code{\link{rStrauss.sp2}},
  \code{\link{rThomas.sp2}} 
}
\examples{
rM2 <- rMaternII.sp2(kappa=20, r=pi/10, win=sphwin())
rM2
}
\keyword{spatial}
\keyword{datagen}

