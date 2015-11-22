\name{methods.sphppm}
\alias{methods.sphppm}
\alias{anova.sphppm}
\alias{coef.sphppm}
\alias{fitted.sphppm}
\alias{print.sphppm}
\title{
  Methods for Fitted Models on a Sphere
}
\description{
  Methods for standard operations
  for the class of point process models on a sphere.
}
\usage{
 \method{anova}{sphppm}(object, \dots)

 \method{coef}{sphppm}(object, \dots)

 \method{fitted}{sphppm}(object, \dots)

 \method{print}{sphppm}(x, \dots)
}
\arguments{
  \item{object,x}{
    A fitted model of class \code{"sphppm"}
    obtained from \code{\link{sphppm}}
  }
  \item{\dots}{
    Additional arguments passed to other methods.
  }
}
\details{
  These are methods for the generic functions
  \code{\link[stats]{anova}},
  \code{\link[stats]{coef}},
  \code{\link[stats]{fitted}}
  and 
  \code{\link[base]{print}},
  for the class \code{"sphppm"} of fitted point process models
  on a sphere.
}
\value{
  Same as for the generics.
}
\author{
  Adrian Baddeley
}
\seealso{
  \code{\link{sphppm}}
}
\examples{
   W <- sphwin(type="band", param=c(pi/3, 2*pi/3), ref=c(0,0))
   Z <- runif.sphwin(20, W)
   model <- sphppm(Z ~ cos(theta) + sin(phi))
  model
  fitted(model)
  anova(model)

   # likelihood ratio test
   model0 <- sphppm(Z ~ 1)
   anova(model0, model, test="Chi")
}
\keyword{spatial}
\keyword{models}
\keyword{methods}
