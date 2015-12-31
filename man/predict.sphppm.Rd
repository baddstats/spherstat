\name{predict.sphppm}
\alias{predict.sphppm}
\title{
  Compute Predictions of Fitted Point Process Model on Sphere
}
\description{
  Given a point process model fitted to spherical point pattern data,
  this function computes predictions of the model such as the fitted
  intensity.
}
\usage{
\method{predict}{sphppm}(object,
                         newdata = NULL,
                         type = c("intensity", "link", "terms"), \dots)
}
\arguments{
  \item{object}{
    The fitted model. An object of class \code{"sphppm"}
    usually obtained from \code{\link{sphppm}}.
  }
  \item{newdata}{
    Optional. New data specifying locations where predictions should be
    computed. Either a spherical point pattern (class \code{"sp2"} or
    \code{"sp3"}) or a data frame. See Details.
  }
  \item{type}{
    Type of prediction: \code{"intensity"} to compute the fitted
    intensity, \code{"link"} to compute the fitted log intensity,
    or \code{"terms"} to evaluate each of the covariate terms in the
    model. 
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is the method for the generic function
  \code{\link[stats]{predict}}
  for the class \code{"sphppm"} of fitted point process models on a
  sphere.

  Without additional arguments, \code{predict(object)} is
  equivalent to \code{fitted(object)} and computes the fitted intensity
  at each of the data points to which the model was originally fitted.

  By specifying \code{newdata} the user can evaluate the predictions of
  the model at other locations. Here \code{newdata} should be either a
  spherical point pattern (class \code{"sp2"} or \code{"sp3"}) or a data frame
  with columns named \code{"theta"}, \code{"phi"}, \code{"x1"},
  \code{"x2"} and \code{"x3"} giving the
  angular and Cartesian coordinates of the locations where prediction is
  required. 
}
\value{
  Numeric vector or matrix.
}
\author{
  Adrian Baddeley
}
\seealso{
  \code{\link{sphppm}},
  \code{\link{fitted.sphppm}}.
}
\examples{
   W <- sphwin(type="band", param=c(pi/3, 2*pi/3), ref=c(0,0))
   Z <- runif.sphwin(100, W)
   fut <- sphppm(Z ~ x1 + x2 + x3)
   a <- predict(fut)
   a[1:10]
   b <- predict(fut, type="terms")
   b[1:5,]

   U <- runif.sphwin(3, W)
   predict(fut, newdata=U)
   predict(fut, newdata=U, type="terms")
}
\keyword{spatial}
\keyword{models}
\keyword{methods}
