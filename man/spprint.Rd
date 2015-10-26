\name{spprint}
\alias{print.sp2}
\alias{print.sp3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print Brief Details of a Point Pattern Dataset observed on a (Subset of) a Sphere
}
\description{
Prints a very brief description of a point pattern dataset.
}
\usage{
print.sp2(x, ...)
print.sp3(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Point pattern (an object of class \code{sp2} for \code{print.sp2}; an object of class \code{sp3} for \code{print.sp3})
}

\item{...}{
Ignored
}

}
\details{
A very brief description of the point pattern \code{x} is printed. 

This is a method for the generic function \code{\link{print}}. 

}

\author{
Tom Lawrence
}

\note{
These functions are the analogues for point processes on the sphere of the function \code{\link{print.ppp}} in \code{\link{spatstat}}, which is the corresponding function for point processes in R^2.  Hence elements of this help page have been taken from \code{\link{print.ppp}} with the permission of A. J. Baddeley.  This enables the information on this help page to be consistent with that for \code{\link{print.ppp}}.  It is hoped that this will minimise or remove any confusion for users of both \code{\link{spatstat}} and \code{\link{spherstat}}.
}

\seealso{
\code{\link{summary.sp2}}, \code{\link{summary.sp3}}, \code{\link{print.sphwin}}, \code{\link{summary.sphwin}}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{"spatial"}
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
