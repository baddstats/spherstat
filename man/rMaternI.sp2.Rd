\name{rMaternI.sp2}
\alias{rMaternI.sp2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Matern Model I
}
\description{
Generate a random point pattern, a simulated realisation of the Matern Model I inhibition process model. 
}
\usage{
rMaternI.sp2(kappa, r, win = sphwin(type="sphere"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kappa}{
Intensity of the Poisson process of proposal points (a positive number).
}
  \item{r}{
Inhibition distance
}
   \item{win}{
The window in which to simulate the pattern, an object of type \code{\link{sphwin}}.
}
}
\details{
This algorithm generates a realisation of Matern's Model I inhibition process inside the window \code{win}.  First, a simulated stationary Poisson process with intensity \code{kappa} is generated.  The number of proposal points will equal the number of points in this pattern; the original pattern is ignored after this step.

The simulated pattern is comprised of proposal points.  The first proposal point is automatically retained. The remaining proposal points are tested one at a time, and only accepted if the minimum great circle distance between it and the pattern of alreadly accepted proposal points is greater than great circle distance \code{r}.

The final pattern of retained points constitute Matern's Model I.
}
\value{
The simulated point pattern (an object of class \code{\link{sp2}}).
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tom Lawrence tjlawrence@bigpond.com
}
\note{
This function is the analogue for point processes on the sphere of the function \code{\link{rMaternI}} in \code{\link{spatstat}}.  Hence elements of this help page have been taken from that for \code{\link{rMaternI}}, with the permission of A. J. Baddeley.  This enables the information on this help page to be consistent with that for\code{\link{rMaternI}}.  It is hoped that this will minimise or remove any confusion for users of both \code{\link{spatstat}} and \code{\link{spherstat}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{link{rHardcore}}, \code{\link{rMatClust.sp2}}, \code{link{rMaternI}}, \code{\link{rMaternII.sp2}}, \code{\link{rpoispp.sp2}}, \code{\link{rStrauss.sp2}}, \code{\link{rThomas.sp2}}
}
\examples{
rM1 <- rMaternI.sp2(kappa=20, r=pi/10, win=sphwin())
rM1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
