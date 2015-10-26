\name{Hsphere}
\alias{Hsphere}
\alias{Hstab.sph}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Empty Space Function H 
}
\description{
Estimates the empty space function \code{H}(r) from a point pattern in a window of on a sphere, or the variance-stabilised form of this function. If the window is the entire sphere, one of the methods used gives the exact value of the \code{H} function.
}
\usage{
Hsphere(X, refpoints, win = sphwin(type = "sphere"), rvals = seq(from = 0, to = pi, length = 512), ...)
Hstab.sph(X, refpoints, win = sphwin(type = "sphere"), rvals = seq(from = 0, to = pi, length = 512), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The observed point pattern, from which an estimate of \code{H(r)} will be computed.  An object of class \code{\link{sp2}} or \code{\link{sp3}}, or a 2 or 3 column matrix giving the locations of points in spherical coordinates.
}
   \item{refpoints}{
The reference points in the window used to estimate the \code{H} function.  An object of class \code{\link{sp2}} or \code{\link{sp3}}, or a 2 or 3 column matrix giving the locations of points in spherical coordinates.
}  
\item{win}{
The window; an object of type \code{sphwin}.  Only required if \code{X} is an object of class \code{matrix}.
}
  \item{rvals}{
Optional.  Numeric vector. The values of the argument \code{r} at which \code{H(r)} should be evaluated. There is a sensible default, which works well the window is either the sphere or contains a significant proportion of the sphere. First-time users are strongly advised not to specify this argument. See below for important conditions on \code{rvals}.
}

  \item{\dots}{
Optional.  Extra detail that can be passed to \code{\link{eroded.areas.sphwin}}.
}
}
\details{
The spherical contact distribution function of a stationary random set \code{X} is the cumulative distribution function H of the distance from a fixed point in space to the nearest point of \code{X}, given that the point lies outside \code{X}. That is, H(r) equals the probability that \code{X} lies closer than \code{r} units away from the fixed point \code{x}, given that \code{X} does not cover \code{x}. 

Let \code{D = d(x,X)} be the shortest distance from an arbitrary point \code{x} to the set \code{X}. Then the spherical contact distribution function is 

\code{H(r) = P(D <= r | D > 0)}

For a point process, the spherical contact distribution function is the same as the empty space function \code{F} discussed in \code{Fsphere}. 

The argument \code{X} may be a point pattern (object of class "sp2" or "sp3") or a window (object of class "sphwin"). It is assumed to be a realisation of a stationary random set. 

The algorithm first calls distmap to compute the distance transform of X, then computes the Kaplan-Meier and reduced-sample estimates of the cumulative distribution following Hansen et al (1999). If \code{conditional=TRUE} (the default) the algorithm returns an estimate of the spherical contact function \code{H(r)} as defined above. If \code{conditional=FALSE}, it instead returns an estimate of the cumulative distribution function \code{H*(r) = P(D <= r)} which includes a jump at \code{r=0} if \code{X} has nonzero area. 


The variance-stabilising transform of \code{H} is \code{asin(sqrt(H(r)))}.  This is implemented in \code{spherstat} as \code{Hstab.sph}, and is a superior statistic to \code{H} used for inference based on a Monte Carlo or envelope-based test.
}
\value{
An object of class \code{\link{fv}}, see \code{\link{fv.object}}, which can be plotted directly using \code{\link{plot.fv}}.

Essentially a data frame containing some or all of the following columns: 
\item{r}{the values of the argument r at which the function \code{H(r)} has been estimated}
 
\item{rs}{the "reduced sample" (border correction) estimator of \code{H(r)}}
 
\item{km}{the spatial Kaplan-Meier estimator of \code{H(r)}}
 
\item{hazard}{the hazard rate \code{lambda(r)} of \code{H(r)} by the spatial Kaplan-Meier method}
 
\item{raw}{the uncorrected estimate of \code{H(r)}, i.e. the empirical distribution of the distance from a fixed point in the window to the nearest point of \code{X}}
 
\item{han}{the Hanisch correction estimator of \code{H(r)}}
 
\item{theo}{the theoretical value of \code{H(r)} for a stationary Poisson process of the same estimated intensity.}
}

\section{warnings }{The reduced sample (border method) estimator of \code{H} is pointwise approximately unbiased, but need not be a valid distribution function; it may not be a nondecreasing function of \code{r}. Its range is always within \code{[0,1]}. 

The spatial Kaplan-Meier estimator of \code{H} is always nondecreasing but its maximum value may be less than 1. 

The estimate of \code{lambda(r)} returned by the algorithm is an approximately unbiased estimate for the integral of \code{lambda()} over the corresponding histogram cell. It may exhibit oscillations due to discretisation effects. We recommend modest smoothing, such as kernel smoothing with kernel width equal to the width of a histogram cell. 
 


}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tom Lawrence tjlawrence@bigpond.com
}
\note{
This function is the analogue for point processes on the sphere of the function \code{\link{Hest}} in \code{\link{spatstat}}, which is the corresponding function for point processes in R^2.  Hence elements of the code for \code{Hsphere} and help page have been taken from \code{\link{Hest}} with the permission of A. J. Baddeley.  This enables the code to be highly efficient and give corresponding output to, and for the information on this help page to be consistent with that for the function \code{\link{Hest}}.  It is hoped that this will minimise or remove any confusion for users of both \code{\link{spatstat}} and \code{\link{spherstat}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
sph <- sphwin(type="sphere")
sph.pp <- rpoispp.sp2(win=sph, lambda=10)
sph.ref <- rpoispp.sp2(win=sph, lambda=150)
sph.Hest <- Hsphere(X=sph.pp, refpoints=sph.ref)
sph.Hstab <- Hstab.sph(X=sph.pp, refpoints=sph.ref)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
