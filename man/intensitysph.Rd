\name{intensitysph}
\alias{intensitysph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Homogeneous Intensity of a Dataset
}
\description{
A function for computing the intensity of a spatial dataset or spatial point process model, under the assumption of homogeneity.
}
\usage{
intensitysph(X, win, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The locations of the points in the point pattern.  An object of class \code{\link{sp2}} or \code{\link{sp3}}, or a 2 column matrix with locations given in spherical coordinates, or a 3 column matrix with locations given in Cartesian coordinates.
}
  \item{win}{
The window.  An object of type \code{\link{sphwin}}.  Only required if \code{X} is not an object of class \code{\link{sp2}} or \code{\link{sp3}}.
}
\item{\dots}{
Extra information.  Not used now, but probably needed when I expand this to other point processes
}
}
\details{
Estimates the intensity by dividing the number of points in \code{X} by the area of \code{win}.
}
\value{
A numeric.
}

\author{
Tom Lawrence
}
\note{
This function is the analogue for point processes on the sphere of the function \code{\link{intensity}} in \code{\link{spatstat}}, which is the corresponding function for point processes in R^2.  Hence elements of this help page have been taken from \code{\link{intensity}} with the permission of A. J. Baddeley.  This enables the information on this help page to be consistent with that for \code{\link{intensity}}.  It is hoped that this will minimise or remove any confusion for users of both \code{\link{spatstat}} and \code{\link{spherstat}}.
}

\seealso{
\code{\link{rpoispp.sphere}}
}
\examples{
sph <- sphwin(type="sphere")
rps <- rpoispp.sp2(win=sph, lambda=10)
intensitysph(X=rps, win=sph)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
