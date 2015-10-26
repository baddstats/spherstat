\name{area.sphwin}
\alias{area.sphwin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the area on the surface of a sphere
}
\description{
Calculates the area of a sphere or a subset of the sphere (cap, band, wedge or polygon)
}
\usage{
area.sphwin(w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{w}{
The window for which the area is calculated.  An object of type \code{\link{sphwin}}.
}
}
\value{
The area of the sphere is returned as an object of type \code{numeric}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tom Lawrence tjlawrence@bigpond.com
}

\note{
This function is the analogue for point processes on the sphere of the function \code{\link{area.owin}} in \code{\link{spatstat}}, which is the corresponding function for point processes in R^2.
}

\seealso{
\code{\link{eroded.areas.sphwin}}, \code{\link{area.owin}}
}
\examples{
sph <- sphwin(type="sphere")
area.sphwin(w=sph)

sphcap <- sphwin(type="sphere", param=c(0, pi/3))
area.sphwin(w=sphcap)
}