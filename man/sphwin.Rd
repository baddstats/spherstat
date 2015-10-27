\name{sphwin}
\alias{sphwin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a Window
}
\description{
Creates an object of class "sphwin" representing an observation window on a sphere 
}
\usage{
sphwin(type, param, ref = c(0, 0), ref2=NULL, rad=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
The type of window that is to be created.  Must be one of \code{"sphere"}, \code{"band"}, \code{"wedge"}, \code{"polygon"} or \code{"quadrangle"}.  See Details for more information.
}
  \item{param}{
(Required for all windows types except \code{"sphere"}) Parameters used to define the window.  See Details for more information.
}
  \item{ref}{
(Required for all windows types except \code{"sphere"}) Reference information used in defining the window.  See Details for more information.
}

  \item{ref2}{
(Only used for windows of type \code{"polygon"}) Reference information used in defining the window.  See Details for more information.
}

\item{rad}{
The radius of the sphere, by default set to 1.  Must be a positive number.}
}
\details{
An object of class \code{sphwin} is required by many of the other functions in this package to describe the observation window. Although it is required to mathematically define the window, the radius of the sphere is not required here, although it is a separate input required in most functions that require an object of class \code{sphwin}, and for some other spherical geometric and trigonometric functions.

If the window is the entire sphere, use \code{type="sphere"}.  The inputs \code{param} and \code{ref} are not relevant and are not required to define this type of window.

If the window can be defined as that space between two circles of colatitude (or is congruent to such a shape) i.e. is a spherical cap or band, use \code{type="band"}.  \code{param} must then be a length 2 numeric containing the colatitudes of the (congruent) circles of colatitude that border the window, in spherical coordinates (radians), lowest first.  \code{ref} must be a length 2 numeric giving the spherical coordinates of one of the two on the sphere that the vector normal to the two circles of colatitude and passing through the centre of the sphere passes through.

If the window is the complement of a band (i.e. two caps with parallel bases) use \code{bandcomp}.  \code{param} is a length 2 numeric containing in ascending order the circles of colatitude that are congruent to the bases of the caps and \code{ref} is the centre of the upper cap (therefore the point diametrically opposite the centre of the lower cap). 

If the window can be defined as that space between two semicircles of longitude (or is congruent to such a shape) i.e. is a spherical wedge, use \code{type="wedge"}.  \code{param} must be a length 2 numeric, the first value is the dihedral angle of the wedge, the second value is the angular correction (that is, after where appropriate the sphere is rotated such that the vertices of the wedge are \code{c(0,0)} and \code{c(pi,0)}, the longitude of the left hand bounding semi-circle of the wedge).  \code{ref} must be a length 2 numeric containing the spherical coordinates of a vertex of the wedge (and, when appropriate, the vertex rotated to the north pole for the purposes of finding the angular correction).

If the window is a polygon, where all boundary arcs are arcs of great circles, use \code{type=polygon}.  \code{param} must be a 2 column matrix, containing the spherical coordinates of all vertices in the order of traversal and such that the coordinates of the first vertex are given in the first and last rows of the matrix.  \code{ref} should be a numeric with the same length as there are vertices in the polygon; the ith value of this numeric takes value 1 if the ith vertex in \code{param} is a reflex angle and 0 otherwise.  If two consecutive vertices are diametrically opposite, then since they are connected by infinitely many different great circle arcs, \code{sphwin} will not be able to define the window.  An error message will appear asking that a vertex be added between the diametrically opposite vertices so that the window can be defined.  An example of this is included below.  Note also that \code{sphwin} assumes all boundary arcs are no longer than a length of half the circumference of a great circle on the sphere (i.e. \code{<= pi*rad}); if this is not the case for your polygon, add an additional vertex along the "oversized" arcs to create two coplanar ones that are both less than or equal to \code{pi*rad}.  The argument \code{ref2} is a numeric of the same length as \code{ref}, whwre the kth value is the circle of colatitude for which the arc between the kth and (k+1)th vertices is congruent to an arc of.  If they are \code{NULL}, \code{ref} and \code{ref2} are by default set to be numerics that respectively contain only 0s and only pi/2s (respectively denoting that the polygon contains no reflex angles, and that all boundary arcs are arcs of great circles).

If the window can be defined as that space between two circles of colatitude and two semicircles of longitude (in other words, the intersection of a spherical band and a spherical wedge), or is congruent to such a shape, use \code{type=="quadrangle"}.  \code{param} must be a length 4 numeric; the first and second values have the same meaning as \code{param} for \code{type="band"}, the third and fourth values have the same meaning as \code{param} for \code{type=="wedge"}.  By thinking of the window as the intersection of a spherical band and spherical wedge, the reference point \code{ref} is the location of one vertex of that wedge (which is also one end of the axis that is orthogonal to the bounding circles of the spherical band).

No other shapes are supported as yet but may be supported in later versions of this package.
}
\value{
A one sentence description of the window.
}

\author{
Tom Lawrence
}


\seealso{
\code{\link{spherstat-package}} for detail on the spherical coordinate system used in this package
}
\examples{
sph <- sphwin(type="sphere")
sph

sphcap <- sphwin(type="band", param=c(0, pi/3), ref=c(0,0))
sphcap

sphband <- sphwin(type="band", param=c(pi/3, 2*pi/3), ref=c(0,0))
sphband

sphwedge <- sphwin(type="wedge", param=c(pi/3, 0), ref=c(0,0))
sphwedge


# The following call to sphwin at the end of the sequence of code below generates an error because the second and third
# vertices are diametrically opposite, so there are infinitely many arcs between them and hence insufficient information 
# to determinme which one is a boundary of the window

x <- matrix(nrow=4, ncol=2)
x[1,] <- c(0,0)
x[2,] <- c(pi/2, 0)
x[3,] <- c(pi/2, pi)
x[4,] <- c(0,0)
# sphpoly <- sphwin(type="polygon", param=x, ref=rep(0, 3))

# The way to overcome this is to use code like that below, to add an additional vertex between the diametrically opposite
# vertices
x <- rbind(x[1:2,], c(pi/2, pi/2), x[3:4,])
sphpoly <- sphwin(type="polygon", param=x, ref=rep(0, 4))

}
\keyword{spatial}
\keyword{datagen}
