\name{spsummary}
\alias{summary.sp2}
\alias{summary.sp3}
\title{
Summary of a Point Pattern Dataset observed on a (Subset of) a Sphere
}
\description{
Prints a useful summary of a point pattern dataset. 
}
\usage{
\method{summary}{sp2}(object, \dots)

\method{summary}{sp3}(object, \dots)
}
\arguments{
  \item{object}{
    Point pattern (an object of class \code{sp2} for \code{summary.sp2};
    an object of class \code{sp3} for \code{summary.sp3}) 
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
A useful summary of the point pattern \code{object} is printed. 

This is a method for the generic function \code{\link{summary}}. 
}
\author{
Tom Lawrence
}

\note{
These functions are the analogues for point processes on the sphere of
the function \code{\link[spatstat]{summary.ppp}} in \pkg{spatstat},
which is the corresponding function for point processes in \eqn{R^2}{R^2}.  Hence
elements of this help page have been taken from
\code{\link[spatstat]{summary.ppp}} with the permission of
A. J. Baddeley.  This enables the information on this help page to be
consistent with that for \code{\link[spatstat]{summary.ppp}}.  It is
hoped that this will minimise or remove any confusion for users of both
\pkg{spatstat} and \pkg{spherstat}. 
}

\seealso{
\code{\link{print.sp2}}, \code{\link{print.sp3}},
\code{\link{print.sphwin}}, \code{\link{summary.sphwin}} 
}
\examples{
X <- matrix(c(1,0,0, 0,0,1, 0,1,0), nrow=3, ncol=3, byrow=TRUE)
win <- sphwin(type="sphere")
X <- sp3(X=X, win=win)
print(X)
summary(X)

Y <- matrix(c(0,0, pi/2, 0, pi/2, pi, pi, 0), nrow=4, ncol=2, byrow=TRUE)
y <- sp2(X=Y, win=win)
print(Y)
summary(Y)
}
\keyword{spatial}
\keyword{methods}

