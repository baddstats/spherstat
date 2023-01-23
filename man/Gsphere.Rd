\name{Gsphere}
\alias{Gsphere}
\title{
Nearest Neighbour Distance Function G 
}
\description{
Using a number of methods, this function estimates the nearest neighbour distance distribution function \eqn{G(r)} from a point pattern in a window of arbitrary shape on a sphere.  If the window is the entire sphere, one of the methods used gives the exact value of the \eqn{G} function.
}
\usage{
Gsphere(X, win = sphwin(type = "sphere"), r = NULL, ...)

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
    Optional.  Numeric vector. The values of the argument \eqn{r} at
    which \eqn{G(r)} should be evaluated. There is a sensible default,
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
The nearest neighbour distance distribution function (also called the "event-to-event", or the "inter-event" distribution) of a point process \code{X} is the cumulative distribution function \eqn{G} of the distance from a typical random point of \code{X} to the nearest other point of \code{X}. 

An estimate of \eqn{G} derived from a spatial point pattern dataset can be used in exploratory data analysis and formal inference about the pattern (Cressie, 1991; Diggle, 1983; Ripley, 1988). In exploratory analyses, the estimate of \eqn{G} is a useful statistic summarising one aspect of the "clustering" of points. For inferential purposes, the estimate of \eqn{G} is usually compared to the true value of \eqn{G} for a completely random (Poisson) point process, which is 

\deqn{G(r) = 1 - exp( - 2 \pi\rho\lambda (1-\cos(r/\rho))}{G(r) = 1 - exp( - 2 * pi * rho * lambda * (1-cos(r/rho))}

where \eqn{\lambda}{lambda} is the intensity (expected number of points per unit area) and \eqn{\rho}{rho} is the radius of the sphere in which \code{X} is observed. Deviations between the empirical and theoretical \eqn{G} curves may suggest spatial clustering or spatial regularity. 

This algorithm estimates the nearest neighbour distance distribution function \eqn{G} from the point pattern \code{X}. It assumes that \code{X} can be treated as a realisation of a homogeneous random spatial point process on a sphere. The window (subset of the sphere) in which \code{X} is observed may be the any of those that can be defined using \code{\link{sphwin}}.

The estimation of \eqn{G} for a window other than the entire sphere is hampered by edge effects arising from the unobservability of points of the random pattern outside the window. An edge correction is needed to reduce bias (Baddeley, 1998; Ripley, 1988). The edge corrections implemented here are the border method ("reduced sample") estimator, the spatial Kaplan-Meier estimator (Baddeley and Gill, 1997) and the Hanisch estimator (Hanisch, 1984).

The argument \code{r} is the vector of values for the distance \code{r} at which \code{G(r)} should be evaluated. It is also used to determine the breakpoints (in the sense of \code{\link[graphics]{hist}}) for the computation of histograms of distances. The estimators are computed from histogram counts. This introduces a discretisation error which is controlled by the fineness of the breakpoints.  The default value of \code{r} is best suited when the window is a sphere or a significant proportion of a sphere; for smaller subsets of the sphere the default value will still give a valid estimate of the \eqn{G} function but the eroded areas function will reach zero sooner, and the value at which that function reaches 0 is the upper limit of \code{r} for the \eqn{G} function.

First-time users would be strongly advised not to specify \code{r}. However, if it is specified, \code{r} must satisfy \code{r[1] = 0}, have length 512, and \code{max(r)} must be larger than the radius of the largest disc contained in the window. Furthermore, the successive entries of \code{r} must be finely spaced.

The algorithm also returns an estimate of the hazard rate function, \eqn{\lambda(r)}{lambda(r)}, of \eqn{G(r)}. The hazard rate is defined as the derivative 

\deqn{\lambda(r) = - (d/dr) \log(1 - G(r))}{lambda(r) = - (d/dr) log(1 - G(r))}

This estimate should be used with caution as \eqn{G} is not necessarily differentiable. 

The naive empirical distribution of distances from each point of the pattern \code{X} to the nearest other point of the pattern, is the exact value of \eqn{G} when the window is the sphere, and a biased estimate of \eqn{G} otherwise. However, in the latter case it is sometimes useful. It is automatically returned by the algorithm. Care should be taken not to use the uncorrected empirical \eqn{G} as if it were an unbiased estimator of \eqn{G}. 

To simply compute the nearest neighbour distance for each point in the pattern, use \code{\link{nndistsph}}.

If \code{X} is an object of class \code{\link{sp2}} or \code{\link{sp3}}, then its radius is used in this calculation; otherwise, the radius of \code{win} is used.
}
\value{
An object of class \code{fv}, see \code{\link[spatstat.explore]{fv.object}}, which can be plotted directly using \code{\link[spatstat.explore]{plot.fv}}.

Essentially a data frame containing some or all of the following columns: 
\item{r}{the values of the argument r at which the function \eqn{G(r)} has been estimated}
 
\item{rs}{the "reduced sample" (border correction) estimator of \eqn{G(r)}}
 
\item{km}{the spatial Kaplan-Meier estimator of \eqn{G(r)}}
 
\item{hazard}{the hazard rate \eqn{\lambda(r)}{lambda(r)} of \eqn{G(r)} by the spatial Kaplan-Meier method}
 
\item{raw}{the uncorrected estimate of \eqn{G(r)}, i.e. the empirical distribution of the distances from each point in the pattern \code{X} to the nearest other point of the pattern}
 
\item{han}{the Hanisch correction estimator of \eqn{G(r)}}

\item{rs.modif}{the "modified reduced sample" (modified border correction) estimator of \eqn{G(r)}}
 
\item{theo}{the theoretical value of \eqn{G(r)} for a stationary Poisson process of the same estimated intensity.}
}

\section{Warnings}{The function \eqn{G} does not necessarily have a density. Any valid c.d.f. may appear as the nearest neighbour distance distribution function of a stationary point process. 

The reduced sample estimator of \eqn{G} is pointwise approximately unbiased, but need not be a valid distribution function; it may not be a nondecreasing function of \code{r}. Its range is always within \code{[0,1]}. 

The spatial Kaplan-Meier estimator of \eqn{G} is always nondecreasing but its maximum value may be less than 1. 
}

\references{

Baddeley, A.J. Spatial sampling and censoring. In O.E. Barndorff-Nielsen, W.S. Kendall and M.N.M. van Lieshout (eds) \emph{Stochastic Geometry: Likelihood and Computation}. Chapman and Hall, 1998. Chapter 2, pp. 37--78. 

Baddeley, A J. and Gill, R.D. (1997) Kaplan-Meier estimators of interpoint distance distributions for spatial point processes. \emph{Annals of Statistics} \bold{25} pp.263--292. 

Bahcall, J. N. and Soneira, R. M. (1981) The distribution of stars to $V=16$th magnitude near the north galactic pole: normalization, clustering properties and counts in various bands \emph{Astrophysical Journal}, \bold{246}, pp. 122--135.

Cressie, N.A.C. Statistics for spatial data. John Wiley and Sons, 1991. 

Diggle, P.J. Statistical analysis of spatial point patterns. Academic Press, 1983. 

Hanisch, K.-H. (1984) Some remarks on estimators of the distribution function of nearest-neighbour distance in stationary spatial point patterns. \emph{Mathematische Operationsforschung und Statistik, series Statistics} \bold{15}, pp.409--412. 

Lawrence, T.J. (2017) Master's Thesis, University of Western Australia.

Ripley, B.D. \emph{Statistical inference for spatial processes}. Cambridge University Press, 1988. 

Robeson, S., Li, A., and Huang C. (2014) Point pattern analysis on the sphere \emph{Spatial Statistics}, \bold{10}, pp. 76--86.

Scott, D. and Tout C. A. (1989) Nearest neighbour analysis of random distributions on a sphere. \emph{Mon. Not. R. astr. Soc.} \bold{241}, pp. 109--117.

Stoyan, D, Kendall, W.S. and Mecke, J. Stochastic geometry and its applications. 2nd edition. Springer Verlag, 1995.
}

\author{
Tom Lawrence <email:tjlawrence@bigpond.com>
}

\note{
This function is the analogue for point processes on the sphere of the function \code{\link[spatstat.explore]{Gest}} in \pkg{spatstat}, which is the corresponding function for point processes in \eqn{R^2}{R^2}.  Hence elements of the code for \code{\link[spatstat.explore]{Gest}} and this help page have been taken from \code{\link[spatstat.explore]{Gest}} with the permission of A. J. Baddeley.  This enables the code to be highly efficient and give corresponding output to, and for the information on this help page to be consistent with that for the function \code{\link[spatstat.explore]{Gest}}.  It is hoped that this will minimise or remove any confusion for users of both \pkg{spatstat} and \pkg{spherstat}.
}

\seealso{
\code{\link{Fsphere}}, \code{\link{Jsphere}}, \code{\link{Ksphere}}, \code{\link[spatstat.explore]{Gest}}
}

\examples{
sph <- sphwin(type="sphere")
sph.pp <- rpoispp.sphwin(win=sph, lambda=10)
sph.Gest <- Gsphere(X=sph.pp)
}
\keyword{spatial}
\keyword{nonparametric}
