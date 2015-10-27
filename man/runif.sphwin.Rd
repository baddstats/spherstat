\name{runif.sphwin}
\alias{runif.sphwin}
\alias{rpoispp.sp2}
\title{
  Generate the uniform distribution or a homogeneous Poisson Point
  Pattern on the sphere 
}
\description{
  Generates the uniform distribution or a random point pattern using the
  (homogeneous) Poisson process, on a sphere or a subset of a
  sphere. Includes CSR (complete spatial randomness).  
}
\usage{
runif.sphwin(n, win = sphwin(type = "sphere"), as.sp=TRUE, sp.dim="2")

rpoispp.sp2(lambda, win = sphwin(type = "sphere"), lmax=NULL, as.sp=TRUE, sp.dim="2", ...)
}
\arguments{
  \item{n}{
    The number of points to be simulated.  An integer.
  }
  \item{lambda}{
    The intensity of the point process to be simulated.  A numeric of
    length 1 or an object of class \code{function}. 
  }
  \item{win}{
    The window, an object of type \code{sphwin}.
  }
  \item{lmax}{
    A numeric of length 1 giving the maximum value of \code{lambda}, if
    \code{lambda} is a function. 
  }
  \item{as.sp}{
    Logical.  If TRUE, returns an object of class as defined by
    \code{sp.dim}.  Otherwise, returns a matrix.  See Value. 
  }
  \item{sp.dim}{
    A string.  Specifies the object that the function returns, if
    \code{as.sp=TRUE}. 
  }
  \item{\dots}{
    Arguments passed to \code{lambda} if it is a function.
  }
}
\details{
  \code{runif.sphwin} generates a realisation of \code{n} points of the
  uniform distribution inside the window \code{win} (if \code{win} is an
  object of class \code{sphwin}). 

  If \code{lambda} is a single number, then \code{rpoispp.sp2} this
  algorithm generates a realisation of the uniform Poisson process (also
  known as Complete Spatial Randomness, CSR) inside the window
  \code{win} with intensity \code{lambda} (points per unit area).
  
  If \code{lambda} is a function, then \code{rpoispp.sp2} uses the
  adaptation to point patterns on the sphere of Lewis and Shedler's (1979)
  rejection sampling algorithm.  According to this algorithm, a uniform
  Poisson process of intensity \code{lmax} is generated, and a point at
  location \code{s} is retained with probability \code{lambda(s)/lmax}.
}
\value{
If \code{as.sp=FALSE}, a two column matrix giving the spherical coordinates of the simulated data.
If \code{as.sp=TRUE} and \code{sp.dim="2"}, an object of class \code{\link{sp2}} giving the simulated data as a point pattern.
If \code{as.sp=TRUE} and \code{sp.dim="3"}, an object of class \code{\link{sp3}} giving the simulated data as a point pattern.
}
\references{
  Lewis, P.A.W. and Shedler, G.S. (1979)
  Simulation of non-homogeneous Poisson processes by thinning.
  \emph{Naval Logistics Quarterly} \bold{26}, 406--413.
}
\author{
Tom Lawrence tjlawrence@bigpond.com
}
\note{
  The function \code{rpoispp.sp2} is the analogue for point processes on
  the sphere of the function \code{\link[spatstat]{rpoispp}} in
  \pkg{spatstat}, which is the corresponding function for point
  processes in R^2.  Hence elements of this help page have been taken
  from \code{\link[spatstat]{rpoispp}} with the permission of
  A. J. Baddeley.  This enables the relevant information on this help
  page to be consistent with that for \code{\link[spatstat]{rpoispp}}.
  It is hoped that this will minimise or remove any confusion for users
  of both \pkg{spatstat} and \pkg{spherstat}. 
}
\seealso{
  \code{\link[spatstat]{rpoispp}}
}
\examples{
sph <- sphwin(type="sphere")
rus <- runif.sphwin(n=150, win=sph, as.sp=TRUE, sp.dim="2")
rus

rps <- rpoispp.sp2(win=sph, lambda=150, as.sp=TRUE, sp.dim="2")
rps
}
\keyword{spatial}
\keyword{datagen}

