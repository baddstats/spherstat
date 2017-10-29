\name{rMatClust.sphwin}
\alias{rMatClust.sphwin}
\title{
  Simulate the analogue of the Matern Cluster Process on the sphere
}
\description{
  Generate a random point pattern, a simulated realisation of the
  analogue of the Matern Cluster Process on the sphere.  
}
\usage{
  rMatClust.sphwin(kappa, scale, mu, win=sphwin(type="sphere"), parents=FALSE, nsim=1, drop=TRUE, expand=TRUE, as.sp=TRUE, ndim="2")
}
\arguments{
  \item{kappa}{
    Intensity of the Poisson process of cluster centres. A single
    positive number. 
  }
  \item{scale}{
    Radius parameter of the clusters.
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
  \item{nsim}{
    Number of simulated realisations to be generated.
  }
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the result will be a point pattern, rather than a list containing a point pattern. 
  }
  \item{expand}{ 
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
  This algorithm generates a realisation of the analogue on the sphere of
  Matern's cluster process, a special case of the Neyman-Scott process,
  inside the window \code{win}.  

  This algorithm generates a uniform Poisson point process of ``parent''
  points with intensity \code{kappa}. Then each parent point is replaced
  by a random cluster of ``offspring'' points, the number of points per
  cluster being Poisson (\code{mu}) distributed, and their positions
  being placed and uniformly inside a disc of radius \code{r} centred on
  the parent point. The resulting point pattern is a realisation of the
  classical ``stationary Matern cluster process'' generated inside the
  window \code{win}. This point process has intensity \code{kappa * mu}. 
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
Lawrence, T.J. (2017) Master's Thesis, University of Western Australia.

Matern, B. (1960) \emph{Spatial Variation}. Meddelanden fraan Statens Skogsforskningsinstitut, volume 59, number 5. Statens Skogsforskningsinstitut, Sweden. 

Matern, B. (1986) \emph{Spatial Variation}. Lecture Notes in Statistics 36, Springer-Verlag, New York. 

Waagepetersen, R. (2007) An estimating function approach to inference for inhomogeneous Neyman-Scott processes. \emph{Biometrics} \bold{63}, 252--258. 

}
\author{
  Tom Lawrence <email:tjlawrence@bigpond.com>
  and Adrian Baddeley
}
\note{
  This function is the analogue for point processes on the sphere of the
  function \code{\link[spatstat]{rMatClust}} in \pkg{spatstat}, which is the
  corresponding function for point processes in \eqn{R^2}{R^2}.  Hence elements of
  this help page have been taken from that for \pkg{spatstat} with
  the permission of A. J. Baddeley. This enables the information on this
  help page to be consistent with that for
  \code{\link[spatstat]{rMatClust}}.  It is hoped that this will
  minimise or remove any confusion for users of both 
  \pkg{spatstat} and \pkg{spherstat}. 
}
\seealso{
  \code{\link{rHardcore.sphwin}}, \code{\link[spatstat]{rMatClust}},
  \code{\link{rMaternI.sphwin}}, \code{\link{rMaternII.sphwin}},
  \code{\link{rpoispp.sphwin}}, \code{\link{rStrauss.sphwin}},
  \code{\link{rThomas.sphwin}} 
}
\examples{
rM1 <- rMatClust.sphwin(250, 0.04*pi, 5, win=sphwin(), parents=FALSE)
rM1
}
\keyword{spatial}
\keyword{datagen}
