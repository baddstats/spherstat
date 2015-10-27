\name{rThomas.sp2}
\alias{rThomas.sp2}
\title{
  Simulate analogue of the Thomas Process on the sphere
}
\description{
  Generate a random point pattern using the (homogeneous) Poisson
  process, or uniform distribution, on a sphere or a subset of a
  sphere. Includes CSR (complete spatial randomness).  
}
\usage{
rThomas.sp2(kappa, sigma, mu, win, parents=FALSE)
}
\arguments{
  \item{kappa}{
    Intensity of the Poisson process of cluster centres. A single
    positive number. 
  }
  \item{sigma}{
    The intensity of the point process to be simulated.
  }
  \item{mu}{
    Mean number of points per cluster (a single positive number)
  }
  \item{win}{
    Window in which to simulate the pattern.  An object of class \code{sphwin}.
  }
  \item{parents}{
    Logical.  If \code{TRUE}, the parent points are included in the
    output, if \code{FALSE} parents are not included. 
  }
}
\details{
  This algorithm generates a realisation of the analogue on the sphere
  (`modified') Thomas process, a special case of the Neyman-Scott process,
  inside the window \code{win}.  

  In this case, \code{kappa} and \code{mu} are single numbers, and the
  algorithm generates a uniform Poisson point process of ``parent''
  points with intensity \code{kappa}. Then each parent point is replaced
  by a random cluster of ``offspring'' points, the number of points per
  cluster being Poisson (\code{mu}) distributed, and their positions
  being Fisher displacements from the cluster parent location. The
  resulting point pattern is a realisation of the analogue on the sphere
  of the classical ``stationary Thomas process'' generated inside the
  window \code{win}. This point process has intensity \code{kappa * mu}.  
}
\value{
  An object of class \code{"sp2"}, giving the locations of points in
  the simulation. 
}
\references{
  Lawrence, T.J. (2016) Master's Thesis, University of Western Australia.
}
\author{
Tom Lawrence <email:tjlawrence@bigpond.com>
}
\note{
  This function is the analogue for point processes on the sphere of the
  function \code{\link[spatstat]{rThomas}} in \pkg{spatstat}, which is
  the corresponding function for point processes in R^2.  Hence elements
  of this help page have been taken from that for
  \code{\link[spatstat]{rThomas}} with the permission of A. J. Baddeley.
  This enables the code to be highly efficient and give corresponding
  output to \code{\link[spatstat]{rThomas}}, and for the information on
  this help page to be consistent with that for
  \code{\link[spatstat]{rThomas}}.  It is hoped that this will minimise
  or remove any confusion for users of both \pkg{spatstat} and
  \pkg{spherstat}. 
}
\seealso{
  \code{\link{rFisher}}, \code{link{rHardcore.sp2}},
  \code{\link{rMatClust.sp2}}, \code{\link{rMaternI.sp2}},
  \code{\link{rMaternII.sp2}}, \code{\link{rpoispp.sp2}},
  \code{\link{rStrauss.sp2}}, \code{\link[spatstat]{rThomas}} 
}
\examples{
rT1 <- rThomas.sp2(250, 25, 5, sphwin(), FALSE)
rT1
}
\keyword{spatial}
\keyword{datagen}

