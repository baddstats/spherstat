\name{Ksphere}
\alias{Ksphere}
\title{
  Estimate of Spherical K Function
}
\description{
  Given a point pattern on the sphere, or on a region of the sphere,
  compute an estimate of the \eqn{K}-function.
}
\usage{
  Ksphere(X, win=sphwin(), rvals,
          correction=c("un", "iso", "rs", "rsm"),
          ratio=FALSE)
}
\arguments{
  \item{X}{
    The observed point pattern, from which an estimate of \eqn{K(r)} will be
    computed.  An object of class \code{\link{sp2}} or \code{\link{sp3}}, or
    a 2 or 3 column matrix giving the locations of points in spherical
    coordinates. 
  }
  \item{win}{
    The window; an object of type \code{sphwin}.  Only required if
    \code{X} is an object of class \code{matrix}. 
  }
  \item{rvals}{
    Optional.  Numeric vector. The values of the argument \eqn{r} at
    which \code{K(r)} should be evaluated.  
  }
  \item{correction}{
    Character string or character vector specifying the choice of
    edge corrections.
  }
  \item{ratio}{
    Logical value specifying whether to store information about the
    numerator and denominator of the estimator, for use by
    other functions.
  }
}
\details{
  To be written.
}
\value{
  An object of class \code{"fv"}, see \code{\link[spatstat]{fv.object}},
  which can be plotted directly using \code{\link[spatstat]{plot.fv}}. 
}
\author{
  Tom Lawrence.
}
\seealso{
  \code{\link{Fsphere}}
}
\examples{
   win <- sphwin()
   X <- runif.sphwin(100, win)
   rv <- seq(0, pi/2, length=512)
   K <- Ksphere(X, win, rv)
   K
   plot(K)

   sphcap <- sphwin(type="band", param=c(0, pi/3), ref=c(0,0))
   Xcap <- runif.sphwin(100, sphcap)
   rv3 <- seq(0, pi/3, length=512)
   Kcap <- Ksphere(Xcap, rvals=rv3)
   plot(Kcap, xlim=c(0, pi/4))

   sphband <- sphwin(type="band", param=c(pi/3, 2*pi/3), ref=c(0,0))
   Xb <- runif.sphwin(100, sphband)
   Kb <- Ksphere(Xb, rvals=rv3)
   plot(Kb)

   sphbc <- sphwin(type="bandcomp", param=c(pi/3, 2*pi/3), ref=c(0,0))
   Xbc <- runif.sphwin(100, sphbc)
   Kbc <- Ksphere(Xbc, rvals=rv3)
   plot(Kbc)
}
\keyword{spatial}
\keyword{nonparametric}
