\name{simulate.sphppm}
\alias{simulate.sphppm}
\title{
  Simulate Realisations of a Fitted Point Process Model on Sphere
}
\description{
  Given a point process model fitted to spherical point pattern data,
  this function generates one or more simulated realisations of the model.
}
\usage{
\method{simulate}{sphppm}(object, nsim = 1, \dots, win, drop = TRUE)
}
\arguments{
  \item{object}{
    The fitted model. An object of class \code{"sphppm"}
    usually obtained from \code{\link{sphppm}}.
  }
  \item{win}{
    The window in which to generate the simulated realisations.
    An object of class \code{"sphwin"}.
    Default is the window for the data to which the model was originally fitted.
  }
  \item{nsim}{
    Number of simulations. Integer.
  }
  \item{\dots}{
    Ignored.
  }
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE}, the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\details{
  This is the method for the generic function
  \code{\link[stats]{simulate}}
  for the class \code{"sphppm"} of fitted point process models on a
  sphere.

  Simulated realisations of the model are generated using
  \code{\link{rpoispp.sphwin}}, using \code{\link{predict.sphppm}} to
  evaluate the intensity function.
}
\value{
  If \code{nsim=1} and \code{drop=TRUE} (the defaults),
  a spherical point pattern (object of class \code{"sp2"}).
  Otherwise, a list of point patterns.
}
\author{
  Adrian Baddeley.
}
\seealso{
  \code{\link{sphppm}}, \code{\link{rpoispp.sphwin}}
}
\examples{
   W <- sphwin(type="band", param=c(pi/3, 2*pi/3), ref=c(0,0))
   Z <- runif.sphwin(100, W)
   fut <- sphppm(Z ~ x1 + x2 + x3)
   simulate(fut, 2)
}
\keyword{spatial}
\keyword{datagen}

