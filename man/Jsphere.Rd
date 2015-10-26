\name{Jsphere}
\alias{Jsphere}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the J-function 
}
\description{
Estimates the summary function \code{J(r)} for a point pattern in a window of arbitrary shape on a (subset of) a sphere.
}
\usage{
Jsphere(X, refpoints, rvals = seq(from = 0, to = pi, length = 512), ..., eps = NULL, r = rvals, breaks=NULL, correction=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The observed point pattern, from which an estimate of \code{J(r)} will be computed.  An object of class \code{\link{sp2}} or \code{\link{sp3}}, or a 2 or 3 column matrix giving the locations of points in spherical coordinates.
}
   \item{refpoints}{
The reference points in the window used to estimate the F function (a prerequisite to calculating the \code{J(r)}, see Details section).  An object of class \code{\link{sp2}} or \code{\link{sp3}}, or a 2 or 3 column matrix giving the locations of points in spherical coordinates.
}  
  \item{rvals}{
Optional.  Numeric vector. The values of the argument \code{r} at which \code{F(r)} should be evaluated. There is a sensible default, which works well the window is either the sphere or contains a significant proportion of the sphere. First-time users are strongly advised not to specify this argument. See below for important conditions on \code{rvals}.
}

  \item{\dots}{
Optional.  Extra detail that can be passed to \code{\link{eroded.areas.sphwin}}.
}

\item{eps}{
Internal argument
}

\item{r}{
Internal argument
}

\item{breaks}{
Internal argument
}

\item{correction}{
The corrections to be applied.  If set to \code{NULL}, all estimators are given.  Specific estimators can be called using the arguments \code{un}, \code{raw}, \code{rs}, \code{rs.modif}, \code{cs}, \code{han} and \code{km}.
}

}
\details{
The \code{J} function (Van Lieshout and Baddeley, 1996) of a stationary point process is defined as 

 \code{J(r) = (1-G(r))/(1-F(r))}

where \code{G(r)} is the nearest neighbour distance distribution function of the point process (see \code{Gsphere}) and \code{F(r)} is its empty space function (see \code{Fsphere}). 

For a completely random (uniform Poisson) point process, the \code{J}-function is identically equal to 1. Deviations \code{J(r) < 1} or \code{J(r) > 1} typically indicate spatial clustering or spatial regularity, respectively. The J-function is one of the few characteristics that can be computed explicitly for a wide range of point processes. See Van Lieshout and Baddeley (1996), Baddeley et al (2000), Thonnes and Van Lieshout (1999) for further information. 

An estimate of \code{J} derived from a spatial point pattern dataset can be used in exploratory data analysis and formal inference about the pattern. The estimate of \code{J(r)} is compared against the constant function 1. Deviations \code{J(r) < 1} or \code{J(r) > 1} may suggest spatial clustering or spatial regularity, respectively. 

This algorithm estimates the \code{J}-function from the point pattern \code{X}. It assumes that \code{X} can be treated as a realisation of a stationary (spatially homogeneous) random spatial point process in the plane, observed through a bounded window. The window (which is specified in \code{X} as \code{X$win}) may have arbitrary shape. 

The argument \code{X} is interpreted as a point pattern object (of class \code{sp2} or \code{sp3}). 

The functions \code{Fsphere} and \code{Gsphere} are called to compute estimates of \code{F(r)} and \code{G(r)} respectively. These estimates are then combined by simply taking the ratio \code{J(r) = (1-G(r))/(1-F(r))}. 

In fact several different estimates are computed using different edge corrections (Baddeley, 1998). 

The Kaplan-Meier estimate (returned as \code{km}) is the ratio \code{J = (1-G)/(1-F)} of the Kaplan-Meier estimates of \code{1-F} and \code{1-G} computed by \code{Fsphere} and \code{Gsphere} respectively. This is computed if \code{correction=NULL} or if correction includes "\code{km}". 

The Hanisch-style estimate (returned as \code{han}) is the ratio \code{J = (1-G)/(1-F)} where \code{F} is the Chiu-Stoyan estimate of \code{F} and \code{G} is the Hanisch estimate of \code{G}. This is computed if \code{correction=NULL} or if correction includes "\code{cs}" or "\code{han}". 

The reduced-sample or border corrected estimate (returned as \code{rs}) is the same ratio \code{J = (1-G)/(1-F)} of the border corrected estimates. This is computed if \code{correction=NULL} or if correction includes "\code{rs}" or "\code{border}". 

These edge-corrected estimators are slightly biased for \code{J}, since they are ratios of approximately unbiased estimators. The logarithm of the Kaplan-Meier estimate is exactly unbiased for \code{log J}. 

The uncorrected estimate (returned as \code{un} and computed only if correction includes "\code{none}") is the ratio \code{J = (1-G)/(1-F)} of the uncorrected ("raw") estimates of the survival functions of\code{F} and \code{G}, which are the empirical distribution functions of the empty space distances \code{Fsphere(X,...)$raw} and of the nearest neighbour distances \code{Gsphere(X,...)$raw}. The uncorrected estimates of \code{F} and \code{G} are severely biased. However the uncorrected estimate of \code{J} is approximately unbiased (if the process is close to Poisson); it is insensitive to edge effects, and should be used when edge effects are severe (see Baddeley et al, 2000). 

The algorithm for \code{Fsphere} uses two discrete approximations which are controlled by the parameter \code{eps} and by the spacing of values of \code{r} respectively. See \code{Fsphere} for details. First-time users are strongly advised not to specify these arguments. 

Note that the value returned by \code{Jsphere} includes the output of \code{Fsphere} and \code{Gsphere} as attributes (see the last example below). If the user is intending to compute the \code{F}, \code{G} and \code{J} functions for the point pattern, it is only necessary to call \code{Jsphere}.

}
\value{
An object of class \code{\link{fv}}, see \code{\link{fv.object}}, which can be plotted directly using \code{\link{plot.fv}}.

Essentially a data frame containing some or all of the following columns: 
\item{r}{the values of the argument r at which the function \code{J(r)} has been estimated}
 
\item{rs}{the "reduced sample" (border correction) estimator of \code{J(r)}}
 
\item{km}{the spatial Kaplan-Meier estimator of \code{J(r)}}
 
\item{hazard}{the hazard rate \code{lambda(r)} of \code{J(r)} by the spatial Kaplan-Meier method}
 
\item{raw}{the uncorrected estimate of \code{J(r)}, i.e. the empirical distribution of the distance from a fixed point in the window to the nearest point of \code{X}}
 
\item{han}{the Hanisch correction estimator of \code{J(r)}}
 
\item{theo}{the theoretical value of \code{J(r)} for a stationary Poisson process of the same estimated intensity.}
The data frame also has attributes:
\item{F}{the output of \code{Fsphere} for this point pattern, containing three estimates of the empty space function \code{F(r)} and an estimate of its hazard function} 
 
\item{G}{
the output of \code{Gsphere} for this point pattern, containing three estimates of the nearest neighbour distance distribution function \code{G(r)} and an estimate of its hazard function}
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
This function is the analogue for point processes on the sphere of the function \code{\link{Jest}} in \code{\link{spatstat}}, which is the corresponding function for point processes in R^2.  Hence elements of the code for \code{Jsphere} and help page have been taken from \code{\link{Jest}} with the permission of A. J. Baddeley.  This enables the code to be highly efficient and give corresponding output to, and for the information on this help page to be consistent with that for the function \code{\link{Jest}}.  It is hoped that this will minimise or remove any confusion for users of both \code{\link{spatstat}} and \code{\link{spherstat}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
sph <- sphwin(type="sphere")
sph.pp <- rpoispp.sp2(win=sph, lambda=10)
sph.ref <- rpoispp.sp2(win=sph, lambda=150)
sph.Jest <- Jsphere(X=sph.pp, refpoints=sph.ref)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
