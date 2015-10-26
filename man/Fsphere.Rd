\name{Fsphere}
\alias{Fsphere}
\alias{Fstab.sph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Empty Space Function F 
}
\description{
Estimates the summary function \code{F(r)} for a point pattern in a window of arbitrary shape on a (subset of) a sphere.
}
\usage{
Fsphere(X, refpoints, win = sphwin(type = "sphere"), rvals = seq(from = 0, to = pi, length = 512), ...)
Fstab.sph(X, refpoints, win = sphwin(type = "sphere"), rvals = seq(from = 0, to = pi, length = 512), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The observed point pattern, from which an estimate of \code{F(r)} will be computed.  An object of class \code{\link{sp2}} or \code{\link{sp3}}, or a 2 or 3 column matrix giving the locations of points in spherical coordinates.
}
   \item{refpoints}{
The reference points in the window used to estimate the F function.  An object of class \code{\link{sp2}} or \code{\link{sp3}}, or a 2 or 3 column matrix giving the locations of points in spherical coordinates.
}  
\item{win}{
The window; an object of type \code{sphwin}.  Only required if \code{X} is an object of class \code{matrix}.
}
  \item{rvals}{
Optional.  Numeric vector. The values of the argument \code{r} at which \code{F(r)} should be evaluated. There is a sensible default, which works well the window is either the sphere or contains a significant proportion of the sphere. First-time users are strongly advised not to specify this argument. See below for important conditions on \code{rvals}.
}

  \item{\dots}{
Optional.  Extra detail that can be passed to \code{\link{eroded.areas.sphwin}}.
}
}
\details{
The empty space function (also called the "spherical contact distribution" or the "point-to-nearest-event" distribution) of a stationary point process \code{X} is the cumulative distribution function \code{F} of the distance from a fixed point in space to the nearest point of \code{X}. 

An estimate of \code{F} derived from a spatial point pattern dataset can be used in exploratory data analysis and formal inference about the pattern (Cressie, 1991; Diggle, 1983; Ripley, 1988). In exploratory analyses, the estimate of \code{F} is a useful statistic summarising the sizes of gaps in the pattern. For inferential purposes, the estimate of \code{F} is usually compared to the true value of \code{F} for a completely random (Poisson) point process, which is 

\code{F(r) = 1 - exp( - 2 * lambda * pi * (1 - cos (r/rad) ) ) }

where \code{lambda} is the intensity (expected number of points per unit area). Deviations between the empirical and theoretical \code{F} curves may suggest spatial clustering or spatial regularity. 

This algorithm estimates the empty space function \code{F} from the point pattern \code{X}. It assumes that \code{X} can be treated as a realisation of a stationary (spatially homogeneous) random spatial point process in the plane, observed through a bounded window. The window (which is specified in \code{X}) may have arbitrary shape. 

The argument \code{X} is interpreted as a point pattern object (of class "ppp", see \code{\link{ppp.object}}) and can be supplied in any of the formats recognised by \code{\link{as.ppp}}. 

The algorithm uses a discrete approximations which is controlled by the spacing of values of \code{r}. (See below for details.) First-time users are strongly advised not to specify these arguments. 

The estimation of \code{F} is hampered by edge effects arising from the unobservability of points of the random pattern outside the window. An edge correction is needed to reduce bias (Baddeley, 1998; Ripley, 1988). The edge corrections implemented here are the border method or "reduced sample" estimator, the spatial Kaplan-Meier estimator (Baddeley and Gill, 1997) and the Chiu-Stoyan estimator (Chiu and Stoyan, 1998). 

The argument \code{r} is the vector of values for the distance \code{r} at which \code{F(r)} should be evaluated. It is also used to determine the breakpoints (in the sense of \code{hist} for the computation of histograms of distances). The estimators are computed from histogram counts. This introduces a discretisation error which is controlled by the fineness of the breakpoints. 

First-time users would be strongly advised not to specify \code{r}. However, if it is specified, \code{r} must satisfy \code{r[1] = 0}, and \code{max(r)} must be larger than the radius of the largest disc contained in the window. Furthermore, the spacing of successive \code{r} values must be very fine. 

The algorithm also returns an estimate of the hazard rate function, \code{lambda(r)}, of \code{F(r)}. The hazard rate is defined by 

\code{ lambda(r) = - (d/dr) log(1 - F(r)) }

The hazard rate of \code{F} has been proposed as a useful exploratory statistic (Baddeley and Gill, 1994). The estimate of lambda(r) given here is a discrete approximation to the hazard rate of the Kaplan-Meier estimator of \code{F}. Note that \code{F} is absolutely continuous (for any stationary point process \code{X}), so the hazard function always exists (Baddeley and Gill, 1997). 

The naive empirical distribution of distances from each location in the window to the nearest point of the data pattern, is a biased estimate of \code{F}. However this is also returned by the algorithm (if \code{correction="none"}), as it is sometimes useful in other contexts. Care should be taken not to use the uncorrected empirical \code{F} as if it were an unbiased estimator of \code{F}. 


}
\value{
An object of class \code{\link{fv}}, see \code{\link{fv.object}}, which can be plotted directly using \code{\link{plot.fv}}.

Essentially a data frame containing some or all of the following columns: 
\item{r}{the values of the argument r at which the function \code{F(r)} has been estimated}
 
\item{rs}{the "reduced sample" (border correction) estimator of \code{F(r)}}
 
\item{km}{the spatial Kaplan-Meier estimator of \code{J(r)}}
 
\item{hazard}{the hazard rate \code{lambda(r)} of \code{F(r)} by the spatial Kaplan-Meier method}
 
\item{raw}{the uncorrected estimate of \code{F(r)}, i.e. the empirical distribution of the distance from a fixed point in the window to the nearest point of \code{X}}
 
\item{han}{the Hanisch correction estimator of \code{F(r)}}
 
\item{theo}{the theoretical value of \code{F(r)} for a stationary Poisson process of the same estimated intensity.}
}

\section{Note}{
Sizeable amounts of memory may be needed during the calculation
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
Tom Lawrence tjlawrence@bigpond.com
}
\note{
This function is the analogue for point processes on the sphere of the function \code{\link{Fest}} in \code{\link{spatstat}}, which is the corresponding function for point processes in R^2.  Hence elements of the code for \code{Fsphere} and help page have been taken from \code{\link{Fest}} with the permission of A. J. Baddeley.  This enables the code to be highly efficient and give corresponding output to, and for the information on this help page to be consistent with that for the function \code{\link{Fest}}.  It is hoped that this will minimise or remove any confusion for users of both \code{\link{spatstat}} and \code{\link{spherstat}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
sph <- sphwin(type="sphere")
sph.pp <- rpoispp.sp2(win=sph, lambda=10)
sph.ref <- rpoispp.sp2(win=sph, lambda=150)
sph.Fest <- Fsphere(X=sph.pp, refpoints=sph.ref)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
