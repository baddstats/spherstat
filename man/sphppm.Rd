\name{sphppm}
\alias{sphppm}
\title{
  Fit Model to Spherical Point Pattern
}
\description{
  Fit a parametric model for the intensity
  to a point pattern on a sphere or subset of the sphere.
}
\usage{
sphppm(formula)
}
\arguments{
  \item{formula}{
    A \code{formula} in the \R language. The left hand side should
    be the name of a point pattern dataset. The right hand side
    specifies the form of the intensity.
  }
}
\details{
  This is analogous to the \pkg{spatstat} function
  \code{\link[spatstat]{ppm}}.

  The argument \code{formula} should be a formula in the \R language.
  The left hand side should be the name of a point pattern dataset,
  which should belong to one of the classes \code{sp2} or \code{sp3}.
  The right hand side of the formula may involve any of the variables
  \code{theta} (colatitude), \code{phi} (longitude)
  and \code{x1, x2, x3} representing three-dimensional spatial location.

  The model is fitted by the logistic regression algorithm of
  Baddeley et al (2014). This procedure uses a random pattern of
  dummy points, so the fitted model coefficients may change slightly
  if the model is re-fitted. 
}
\value{
  An object of class \code{sphppm}.
  See \code{\link{methods.sphppm}} for available methods.
}
\author{
  Adrian Baddeley.
}
\references{
  Baddeley, A. and Coeurjolly, J.-F. and Rubak, E. and Waagepetersen, R. (2014)
  Logistic regression for spatial Gibbs point processes.
  \emph{Biometrika} \bold{101}, 377--392.
}
\seealso{
  \code{\link{methods.sphppm}}
}
\examples{
   W <- sphwin(type="band", param=c(pi/3, 2*pi/3), ref=c(0,0))
   Z <- runif.sphwin(100, W)
   sphppm(Z ~ 1)
   sphppm(Z ~ cos(theta) + cos(phi))
   sphppm(Z ~ abs(x1))
}
\keyword{spatial}
\keyword{models}

