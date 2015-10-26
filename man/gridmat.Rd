\name{gridmat}
\alias{gridmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a grid on the sphere
}
\description{
Generates a grid on the sphere, or a spherical cap, band, wedge, or quadrangle (subset of the sphere bounded by two circles of colatitude and two semcircles of longitude) 
}
\usage{
gridmat(colats, lons, ncolat=100, nlon=100, cellarea=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{colats}{
A length 2 numeric, giving the minimum and maximum colatitudes to be included in the grid (in that order)
}
  \item{lons}{
A length 2 numeric, giving the minimum and maximum longitudes to be included in the grid (in that order)
}
  \item{ncolat}{
Optional.  The number of rows in the grid (each row corresponding to a circle of colatitude).
}
  \item{nlon}{
Optional. The number of columns in the grid (each row corresponding to a semicircle of longitude).
}
\item{cellarea}{
Print the area of the cells in a third column of the output matrix.
}
}
\details{
This function generates a grid on the sphere, or a spherical cap, band, wedge, or quadrangle (subset of the sphere bounded by two circles of colatitude and two semcircles of longitude).  The cells in the grid have equal area, and width (measured in angular distance); the cosine The output is the centres of all the grid cells, which can be used in calculations e.g. to estimate the distribution of the distance from a randomly selected point in a window to the nearst point in a spatial point pattern.  Finer grids (i.e. high values of \code{ncolats, nlons}) heighten the accuracy of such estimates but can take a long time for R to generate those estimates.

Note that the shape of the cells is dependent on the values of the input variables.  Ideally,

For a spherical cap or band, set \code{lons=c(0, 2*pi)} and \code{colats} as appropriate; for a spherical wedge, set \code{colats=c(0, pi)} and \code{lons} as appropriate.  For the entire sphere, set \code{colats=c(0, pi), lons=(0, 2*pi)}.

As the output is in spherical coordinates, the radius of the sphere is not required.  Correctly invoking \code{\link{convert3}} on the result of \code{gridmat} gives the locations of the centres of the cells on the sphere of desired radius, in Cartesian coordinates
The function \code{\link{rot.sphere}} can be used to rotate the grid if required.
}
\value{
A 2 column matrix giving the locations of all points in the grid, in spherical coordinates.  The points are ordered by colatitude, then longitude - both in ascending order.
}

\author{
Tom Lawrence
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
gridsph <- gridmat(colats=c(0, pi), lons=c(0, 2*pi), ncolat=100, nlon=100)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
