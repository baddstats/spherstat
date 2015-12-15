\name{spprint}
\alias{print.sp2}
\alias{print.sp3}
\title{
Print Brief Details of a Point Pattern Dataset observed on a (Subset of) a Sphere
}
\description{
Prints a very brief description of a point pattern dataset.
}
\usage{
\method{print}{sp2}(x, ...)

\method{print}{sp3}(x, ...)
}
\arguments{
  \item{x}{
    Point pattern (an object of class \code{sp2} for \code{print.sp2};
    an object of class \code{sp3} for \code{print.sp3}) 
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
A very brief description of the point pattern \code{x} is printed. 

This is a method for the generic function \code{\link{print}}. 
}
\author{
Tom Lawrence <email: tjlawrence@bigpond.com>
}

\note{
These functions are the analogues for point processes on the sphere of
the function \code{\link[spatstat]{print.ppp}} in \pkg{spatstat}, which
is the corresponding function for point processes in \eqn{R^2}{R^2}.  Hence
elements of this help page have been taken from
\code{\link[spatstat]{print.ppp}} with the permission of A. J. Baddeley.
This enables the information on this help page to be consistent with
that for \code{\link[spatstat]{print.ppp}}.  It is hoped that this will
minimise or remove any confusion for users of both \pkg{spatstat} and
\pkg{spherstat}. 
}

\seealso{
\code{\link{summary.sp2}}, \code{\link{summary.sp3}},
\code{\link{print.sphwin}}, \code{\link{summary.sphwin}} 
}
\examples{
X <- matrix(c(1,0,0, 0,0,1, 0,1,0), nrow=3, ncol=3, byrow=TRUE)
win <- sphwin(type="sphere")
X <- sp3(X=X, win=win)
print(X)

Y <- matrix(c(0,0, pi/2, 0, pi/2, pi, pi, 0), nrow=4, ncol=2, byrow=TRUE)
y <- sp2(X=Y, win=win)
print(Y)
}
\keyword{spatial}
\keyword{methods}
