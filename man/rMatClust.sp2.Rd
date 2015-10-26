\name{rMatClust.sp2}
\alias{rMatClust.sp2}
\alias{rpoispp.sphere}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate the analogue of the Matern Cluster Process on the sphere
}
\description{
Generate a random point pattern, a simulated realisation of the analogue of the Matern Cluster Process on the sphere. 
}
\usage{
rMatClust.sp2(kappa, r, mu, win, parents=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kappa}{
Intensity of the Poisson process of cluster centres. A single positive number.
}
  \item{r}{
Radius parameter of the clusters.
}
  \item{mu}{
Mean number of points per cluster (a single positive number) 
}
  \item{win}{
Window in which to simulate the pattern.  An object of class \code{sphwin}.
}
  \item{parents}{
Logical.  If \code{TRUE}, the parent points are included in the output, if \code{FALSE} parents are not included.
}

}
\details{
This algorithm generates a realisation of the analogue on the sphere of Matern's cluster process, a special case of the Neyman-Scott process, inside the window \code{win}. 

This algorithm generates a uniform Poisson point process of ``parent'' points with intensity \code{kappa}. Then each parent point is replaced by a random cluster of ``offspring'' points, the number of points per cluster being Poisson (\code{mu}) distributed, and their positions being placed and uniformly inside a disc of radius \code{r} centred on the parent point. The resulting point pattern is a realisation of the classical ``stationary Matern cluster process'' generated inside the window win. This point process has intensity \code{kappa * mu}.

}
\value{
An object of class \code{\link{sp2}} giving the locations of the simulated points.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tom Lawrence tjlawrence@bigpond.com
}
\note{

This function is the analogue for point processes on the sphere of the function \code{\link{rMatClust}} in \code{\link{spatstat}}, which is the corresponding function for point processes in R^2.  Hence elements of this help page have been taken from that for \code{\link{spatstat}} with the permission of A. J. Baddeley.  This enables the information on this help page to be consistent with that for \code{\link{rMatClust}}.  It is hoped that this will minimise or remove any confusion for users of both \code{\link{spatstat}} and \code{\link{spherstat}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{link{rHardcore.sp2}}, \code{\link{rMatClust}}, \code{\link{rMaternI.sp2}}, \code{\link{rMaternII.sp2}}, \code{\link{rpoispp.sp2}}, \code{\link{rStrauss.sp2}}, \code{\link{rThomas.sp2}}
}
\examples{
rM1 <- rMatClust(250, 0.04*pi, 5, win=sphwin(), parents=FALSE)
rM1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
