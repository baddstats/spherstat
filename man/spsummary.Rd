\name{spsummary}
\alias{summary.sp2}
\alias{summary.sp3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of a Point Pattern Dataset observed on a (Subset of) a Sphere
}
\description{
Prints a useful summary of a point pattern dataset. 
}
\usage{
summary.sp2(object, ...)
summary.sp3(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Point pattern (an object of class \code{sp2} for \code{summary.sp2}; an object of class \code{sp3} for \code{summary.sp3})
}

\item{...}{
Ignored
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
These functions are the analogues for point processes on the sphere of the function \code{\link{summary.ppp}} in \code{\link{spatstat}}, which is the corresponding function for point processes in R^2.  Hence elements of this help page have been taken from \code{\link{summary.ppp}} with the permission of A. J. Baddeley.  This enables the information on this help page to be consistent with that for \code{\link{summary.ppp}}.  It is hoped that this will minimise or remove any confusion for users of both \code{\link{spatstat}} and \code{\link{spherstat}}.
}

\seealso{
\code{\link{print.sp2}}, \code{\link{print.sp3}}, \code{\link{print.sphwin}}, \code{\link{summary.sphwin}}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{"spatial"}
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
