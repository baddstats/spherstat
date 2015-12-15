\name{rThomas.sphwin}
\alias{rThomas.sphwin}
\title{
  Simulate analogue of the Thomas Process on the sphere
}
\description{
  Generate a random point pattern using the (homogeneous) Poisson
  process, or uniform distribution, on a sphere or a subset of a
  sphere. Includes CSR (complete spatial randomness).  
}
\usage{
rThomas.sphwin(kappa, scale, mu, win=sphwin(type="sphere"), parents=FALSE, nsim=1, drop=TRUE, expand=TRUE, as.sp=TRUE, ndim="2", poisthresh=1e-06)
}
\arguments{
  \item{kappa}{
    Intensity of the Poisson process of cluster centres. A single
    positive number. 
  }
  \item{scale}{
    The concentration parameter of the Fisher distribution that is used the generate the daughter points about each cluster centre.
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
  \item{poisthresh}{
    Numerical threshold below which the model will be treated as a Poisson process. See Details. 
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

If the pair correlation function of the model is very close to that of a Poisson process, deviating by less than \code{poisthresh}, then the model is approximately a Poisson process, and will be simulated as a Poisson process with intensity \eqn{kappa * mu}, using \code{rpoispp.sp}. This avoids computations that would otherwise require huge amounts of memory.  
}
\value{
  If \code{nsim=1} and \code{drop=FALSE} then a single item as described below; otherwise a list containing \code{nsim} items.

  An item is determined by the values of \code{as.sp} and \code{ndim}:
  If \code{as.sp=FALSE} and \code{ndim="2"}, a two column matrix giving the locations of the simulated points.
  If \code{as.sp=FALSE} and \code{ndim="3"}, a three column matrix giving the locations of the simulated points.
  If \code{as.sp=TRUE} and \code{ndim="2"}, an object of class \code{\link{sp2}} giving the locations of the simulated points.
  If \code{as.sp=TRUE} and \code{ndim="3"}, an object of class \code{\link{sp3}} giving the locations of the simulated points.

}
\author{
Tom Lawrence <email:tjlawrence@bigpond.com>
}
\note{
  This function is the analogue for point processes on the sphere of the
  function \code{\link[spatstat]{rThomas}} in \pkg{spatstat}, which is
  the corresponding function for point processes in \eqn{R^2}{R^2}.  Hence elements
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
  \code{\link{rFisher}}, \code{link{rHardcore.sphwin}},
  \code{\link{rMatClust.sphwin}}, \code{\link{rMaternI.sphwin}},
  \code{\link{rMaternII.sphwin}}, \code{\link{rpoispp.sphwin}},
  \code{\link{rStrauss.sphwin}}, \code{\link[spatstat]{rThomas}} 
}
\references{
Diggle, P. J., Besag, J. and Gleaves, J. T. (1976) Statistical analysis of spatial point patterns by means of distance methods. \emph{Biometrics} \bold{32} 659--667. 

Lawrence, T.J. (2016) Master's Thesis, University of Western Australia.

Thomas, M. (1949) A generalisation of Poisson's binomial limit for use in ecology. \emph{Biometrika} \bold{36}, 18--25. 

Waagepetersen, R. (2007) An estimating function approach to inference for inhomogeneous Neyman-Scott processes. \emph{Biometrics} \bold{63}, 252--258. 

}
\examples{
rT1 <- rThomas.sphwin(250, 25, 5, sphwin(), FALSE)
rT1
}
\keyword{spatial}
\keyword{datagen}

