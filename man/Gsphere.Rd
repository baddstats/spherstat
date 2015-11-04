\name{Gsphere}
\alias{Gsphere}
\alias{Gstab.sph}
\title{
Nearest Neighbour Distance Function G 
}
\description{
Using a number of methods, this function estimates the nearest neighbour distance distribution function \code{G(r)} from a point pattern in a window of arbitrary shape on a sphere.  If the window is the entire sphere, one of the methods used gives the exact value of the \code{G} function.
}
\usage{
Gsphere(X, win = sphwin(type = "sphere"), r = NULL, ...)
Gstab.sph(X, win = sphwin(type = "sphere"), r = NULL, ...)
}
\arguments{
  \item{X}{
    The observed point pattern, from which an estimate of \code{G(r)}
    will be computed.  An object of class \code{"sp2"} or
    \code{"sp3"}, or a 2 column matrix giving the locations of points
    in spherical coordinates.
  }
  \item{win}{
    The window; an object of type \code{sphwin}.
}
  \item{r}{
    Optional.  Numeric vector. The values of the argument \code{r} at
    which \code{G(r)} should be evaluated. There is a sensible default,
    which works well the window is either the sphere or contains a
    significant proportion of the sphere. First-time users are strongly
    advised not to specify this argument. 
  }
  \item{\dots}{
    Optional.  Extra detail that can be passed to
    \code{\link{eroded.areas.sphwin}}. 
  }
}
\details{
The nearest neighbour distance distribution function (also called the "event-to-event", or the "inter-event" distribution) of a point process \code{X} is the cumulative distribution function \code{G} of the distance from a typical random point of \code{X} to the nearest other point of \code{X}. 

An estimate of \code{G} derived from a spatial point pattern dataset can be used in exploratory data analysis and formal inference about the pattern (Cressie, 1991; Diggle, 1983; Ripley, 1988). In exploratory analyses, the estimate of \code{G} is a useful statistic summarising one aspect of the "clustering" of points. For inferential purposes, the estimate of \code{G} is usually compared to the true value of G for a completely random (Poisson) point process, which is 

\code{G(r) = 1 - exp( - lambda * pi *2 * (1-cos(r))}

where lambda is the intensity (expected number of points per unit area). Deviations between the empirical and theoretical \code{G} curves may suggest spatial clustering or spatial regularity. 

This algorithm estimates the nearest neighbour distance distribution function \code{G} from the point pattern \code{X}. It assumes that \code{X} can be treated as a realisation of a stationary (spatially homogeneous) random spatial point process in the plane, observed through a bounded window. The window may have arbitrary shape on the sphere. 

The estimation of \code{G} for a window other than the entire sphere is hampered by edge effects arising from the unobservability of points of the random pattern outside the window. An edge correction is needed to reduce bias (Baddeley, 1998; Ripley, 1988). The edge corrections implemented here are the border method ("reduced sample") estimator, the spatial Kaplan-Meier estimator (Baddeley and Gill, 1997) and the Hanisch estimator (Hanisch, 1984).

The argument \code{r} is the vector of values for the distance \code{r} at which \code{G(r)} should be evaluated. It is also used to determine the breakpoints (in the sense of hist) for the computation of histograms of distances. The estimators are computed from histogram counts. This introduces a discretisation error which is controlled by the fineness of the breakpoints.  The default value of \code{r} is best suited when the window is a sphere or a significant proportion of a sphere; for smaller subsets of the sphere the default value will still give a valid estimate of the \code{G} function but the eroded areas function will reach zero sooner, and the value at which that function reaches 0 is the upper limit of r for the \code{G} function.

First-time users would be strongly advised not to specify \code{r}. However, if it is specified, \code{r} must satisfy \code{r[1] = 0}, have length 512, and \code{max(r)} must be larger than the radius of the largest disc contained in the window. Furthermore, the successive entries of \code{r} must be finely spaced. Using the default values of \code{r} for windows other than the entire sphere may impact 

The algorithm also returns an estimate of the hazard rate function, \code{lambda(r)}, of \code{G(r)}. The hazard rate is defined as the derivative 

\code{lambda(r) = - (d/dr) log(1 - G(r))}

This estimate should be used with caution as G is not necessarily differentiable. 

The naive empirical distribution of distances from each point of the pattern \code{X} to the nearest other point of the pattern, is the exact value of \code{G} when the window is the sphere, and a biased estimate of \code{G} otherwise. However, in the latter case it is sometimes useful. It is automatically returned by the algorithm. Care should be taken not to use the uncorrected empirical \code{G} as if it were an unbiased estimator of \code{G}. 

To simply compute the nearest neighbour distance for each point in the pattern, use nndist.sph.

If \code{X} is an object of class \code{\link{sp2}} or \code{\link{sp3}}, then its radius is used in this calculation; otherwise, the radius of \code{win} is used.

The variance-stabilising transform of \code{G} is \code{asin(sqrt(G(r)))}.  This is implemented in \code{spherstat} as \code{Gstab.sph}, and is a superior statistic to \code{G} when used for inference based on a Monte Carlo or envelope-based test.
}
\value{
An object of class \code{\link{fv}}, see \code{\link{fv.object}}, which can be plotted directly using \code{\link{plot.fv}}.

Essentially a data frame containing some or all of the following columns: 
\item{r}{the values of the argument r at which the function \code{G(r)} has been estimated}
 
\item{rs}{the "reduced sample" (border correction) estimator of \code{G(r)}}
 
\item{km}{the spatial Kaplan-Meier estimator of \code{G(r)}}
 
\item{hazard}{the hazard rate \code{lambda(r)} of \code{G(r)} by the spatial Kaplan-Meier method}
 
\item{raw}{the uncorrected estimate of \code{G(r)}, i.e. the empirical distribution of the distances from each point in the pattern X to the nearest other point of the pattern}
 
\item{han}{the Hanisch correction estimator of \code{G(r)}}

\item{rs.modif}{the "modified reduced sample" (modified border correction) estimator of \code{G(r)}}
 
\item{theo}{the theoretical value of \code{G(r)} for a stationary Poisson process of the same estimated intensity.}
}

\section{warnings }{The function \code{G} does not necessarily have a density. Any valid c.d.f. may appear as the nearest neighbour distance distribution function of a stationary point process. 

The reduced sample estimator of \code{G} is pointwise approximately unbiased, but need not be a valid distribution function; it may not be a nondecreasing function of \code{r}. Its range is always within \code{[0,1]}. 

The spatial Kaplan-Meier estimator of \code{G} is always nondecreasing but its maximum value may be less than 1. 


}
\references{
  Lawrence, T.J. (2016) Master's Thesis, University of Western Australia.
}
\author{
Tom Lawrence <email:tjlawrence@bigpond.com>
}
\note{
This function is the analogue for point processes on the sphere of the function \code{\link{Gest}} in \code{\link{spatstat}}, which is the corresponding function for point processes in R^2.  Hence elements of the code for \code{\link{Gest}} and this help page have been taken from \code{\link{Gest}} with the permission of A. J. Baddeley.  This enables the code to be highly efficient and give corresponding output to, and for the information on this help page to be consistent with that for the function \code{\link{Gest}}.  It is hoped that this will minimise or remove any confusion for users of both \code{\link{spatstat}} and \code{\link{spherstat}}.
}
\seealso{
  \code{\link{Fsphere}}
}
\examples{
sph <- sphwin(type="sphere")
sph.pp <- rpoispp.sp2(win=sph, lambda=10)
sph.Gest <- Gsphere(X=sph.pp)
sph.Gstab <- Gstab.sph(X=sph.pp)
}
\keyword{spatial}
\keyword{nonparametric}
