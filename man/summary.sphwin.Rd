\name{summary.sphwin}
\alias{summary.sphwin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of a Spatial Window that is a (Subset of) a Sphere
}
\description{
Prints a useful description of a window object. 
}
\usage{
summary.sphwin(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ 
Window (object of class "\code{sphwin}").
}

\item{...}{
Ignored
}

}
\details{
A useful description of the window \code{object} is printed. 

This is a method for the generic function \code{\link{summary}}. 

}

\author{
Tom Lawrence
}

\note{
This function is the analogue for windows on the sphere of the function \code{\link{summary.owin}} in \code{\link{spatstat}}, which is the corresponding function for point processes in R^2.  Hence elements of this help page have been taken from \code{\link{summary.owin}} with the permission of A. J. Baddeley.  This enables the information on this help page to be consistent with that for \code{\link{summary.owin}}.  It is hoped that this will minimise or remove any confusion for users of both \code{\link{spatstat}} and \code{\link{spherstat}}.
}

\seealso{
\code{\link{summary.sp2}}, \code{\link{summary.sp3}}, \code{\link{print.sp2}}, \code{\link{print.sp2}}, \code{\link{print.sphwin}}
}
\examples{
sph <- sphwin(type="sphere")
print(sph)
summary(sph)

sphcap <- sphwin(type="band", param=c(0, pi/3), ref=c(0,0))
print(sphcap)
summary(sphcap)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{"spatial"}
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
