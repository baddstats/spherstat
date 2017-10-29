\name{Jsphere}
\alias{Jsphere}
\title{
  Estimate the J-function 
}
\description{
  Estimates the summary function \eqn{J(r)} for a point pattern in a
  window of arbitrary shape on a (subset of) a sphere. 
}
\usage{
Jsphere(X, refpoints, r=NULL, \dots, correction=NULL)
}
\arguments{
  \item{X}{
    The observed point pattern, from which an estimate of \eqn{J(r)}
    will be computed.  An object of class \code{\link{sp2}} or
    \code{\link{sp3}}, or a 2 or 3 column matrix giving the locations of
    points in spherical coordinates. 
  }
   \item{refpoints}{
     The reference points in the window used to estimate the \eqn{F} function
     (a prerequisite to calculating \eqn{J(r)}, see Details
     section).  An object of class \code{\link{sp2}} or
     \code{\link{sp3}}, or a 2 or 3 column matrix giving the locations
     of points in spherical coordinates. If \code{NULL}, a homogeneous Poisson process with expected number of points equal to 100 is generated inside the window.
   }  
  \item{r}{
    Optional.  Numeric vector. The values of the argument \eqn{r} at
    which \eqn{F(r)} should be evaluated. There is a sensible default,
    which works well the window is either the sphere or contains a
    significant proportion of the sphere. First-time users are strongly
    advised not to specify this argument.
  }
  \item{\dots}{
    Optional.  Extra detail that can be passed to
    \code{\link{eroded.areas.sphwin}}. 
  }
  \item{correction}{
    The corrections to be applied.  If set to \code{NULL}, all
    estimators are given.  Specific estimators can be called using the
    arguments \code{un}, \code{raw}, \code{rs}, \code{rs.modif},
    \code{cs}, \code{han} and \code{km}. 
  }
}
\details{
  The \eqn{J} function (Van Lieshout and Baddeley, 1996) of a
  stationary point process is defined as  
  \deqn{J(r) = \frac{1-G(r)}{1-F(r)}}{J(r) = (1-G(r))/(1-F(r))}
  where \eqn{G(r)} is the nearest neighbour distance distribution
  function of the point process (see \code{Gsphere}) and \eqn{F(r)} is
  its empty space function (see \code{Fsphere}).  

  For a completely random (uniform Poisson) point process, the
  \eqn{J}-function is identically equal to 1. Deviations \eqn{J(r) < 1} or
  \eqn{J(r) > 1} typically indicate spatial clustering or spatial
  regularity, respectively. The \eqn{J}-function is one of the few
  characteristics that can be computed explicitly for a wide range of
  point processes. See Van Lieshout and Baddeley (1996), Baddeley et al
  (2000), Thonnes and Van Lieshout (1999) for further information.  

  An estimate of \eqn{J} derived from a spatial point pattern dataset
  can be used in exploratory data analysis and formal inference about
  the pattern. The estimate of \eqn{J(r)} is compared against the
  constant function 1. Deviations \eqn{J(r) < 1} or \eqn{J(r) > 1} may
  suggest spatial clustering or spatial regularity, respectively.  

  This algorithm estimates the \eqn{J}-function from the point pattern
  \code{X}. It assumes that \code{X} can be treated as a realisation of a
  stationary (spatially homogeneous) random spatial point process in the
  plane, observed through a bounded window. The window (which is specified
  in \code{X} as \code{X$win}) may have arbitrary shape.  

  The argument \code{X} is interpreted as a point pattern object (of
  class \code{sp2} or \code{sp3}).  

  The functions \code{\link{Fsphere}} and \code{\link{Gsphere}}
  are called to compute
  estimates of \eqn{F(r)} and \eqn{G(r)} respectively. These estimates
  are then combined by simply taking the ratio
  \eqn{J(r) = (1-G(r))/(1-F(r))}.  

  In fact several different estimates are computed using different edge
  corrections (Baddeley, 1998).  

  The Kaplan-Meier estimate (returned as \code{km}) is the ratio
  \eqn{J = (1-G)/(1-F)} of the Kaplan-Meier estimates of \eqn{1-F} and
  \eqn{1-G} computed by \code{Fsphere} and \code{Gsphere}
  respectively. This is computed if \code{correction=NULL} or if
  correction includes \code{"km"}.  

  The Hanisch-style estimate (returned as \code{han}) is the ratio
  \eqn{J = (1-G)/(1-F)} using the Chiu-Stoyan estimate of
  \eqn{F} and the Hanisch estimate of \eqn{G}. This is
  computed if \code{correction=NULL} or if correction includes
  \code{"cs"} or \code{"han"}.  

  The reduced-sample or border corrected estimate (returned as
  \code{rs}) is the same ratio \eqn{J = (1-G)/(1-F)} of the border
  corrected estimates. This is computed if \code{correction=NULL} or if
  correction includes \code{"rs"} or \code{"border"}.  

  These edge-corrected estimators are slightly biased for \eqn{J},
  since they are ratios of approximately unbiased estimators. The
  logarithm of the Kaplan-Meier estimate is exactly unbiased for
  \eqn{log(J)}.  

  The uncorrected estimate (returned as \code{un} and computed only if
  correction includes \code{"none"}) is the ratio
  \eqn{J = (1-G)/(1-F)}
  of the uncorrected ("raw") estimates of the survival functions
  of \eqn{F} and \eqn{G}, which are the empirical distribution
  functions of the empty space distances \code{Fsphere(X,...)$raw} and
  of the nearest neighbour distances \code{Gsphere(X,...)$raw}. The
  uncorrected estimates of \eqn{F} and \eqn{G} are severely
  biased. However the uncorrected estimate of \eqn{J} is approximately
  unbiased (if the process is close to Poisson); it is insensitive to
  edge effects, and should be used when edge effects are severe (see
  Baddeley et al, 2000).  

  The algorithm for \code{Fsphere} uses two discrete approximations
  which are controlled by the parameter \code{eps} and by the spacing of
  values of \code{r} respectively. See \code{Fsphere} for
  details. First-time users are strongly advised not to specify these
  arguments.  

  Note that the value returned by \code{Jsphere} includes the output of
  \code{Fsphere} and \code{Gsphere} as attributes (see the last example
  below). If the user is intending to compute the \eqn{F}, \eqn{G} and
  \eqn{J} functions for the point pattern, it is only necessary to call
  \code{Jsphere}. 
}
\value{
  An object of class \code{"fv"}, see \code{\link[spatstat]{fv.object}},
  which can be plotted directly using \code{\link[spatstat]{plot.fv}}. 

  Essentially a data frame containing some or all of the following columns: 
  \item{r}{the values of the argument r at which the function \eqn{J(r)}
  has been estimated} 
 
  \item{rs}{the "reduced sample" (border correction) estimator of \eqn{J(r)}}
 
  \item{km}{the spatial Kaplan-Meier estimator of \eqn{J(r)}}
 
  \item{hazard}{the hazard rate \eqn{\lambda(r)}{lambda(r)} of \eqn{J(r)} by the
    spatial Kaplan-Meier method} 
 
  \item{raw}{the uncorrected estimate of \eqn{J(r)}, i.e. the empirical
    distribution of the distance from a fixed point in the window to the
    nearest point of \code{X}} 
 
  \item{han}{the Hanisch correction estimator of \eqn{J(r)}}
 
  \item{theo}{the theoretical value of \eqn{J(r)} for a stationary
    Poisson process of the same estimated intensity.} 

  The data frame also has attributes:
  \item{F}{the output of \code{Fsphere} for this point pattern,
    containing three estimates of the empty space function \eqn{F(r)} and
    an estimate of its hazard function}  
 
  \item{G}{
    the output of \code{Gsphere} for this point pattern, containing
    three estimates of the nearest neighbour distance distribution
    function \eqn{G(r)} and an estimate of its hazard function} 
}
\section{Note}{
Sizeable amounts of memory may be needed during the calculation
}
\author{
Tom Lawrence <email:tjlawrence@bigpond.com>
}
\note{
  This function is the analogue for point processes on the sphere of the
  function \code{\link[spatstat]{Jest}} in \pkg{spatstat}, which is the
  corresponding function for point processes in \eqn{R^2}{R^2}.  Hence
  elements of the code for \code{Jsphere} and help page have been taken
  from \code{\link[spatstat]{Jest}} with the permission of
  A. J. Baddeley.  This enables the code to be highly efficient and give
  corresponding output to, and for the information on this help page to
  be consistent with that for the function \code{\link[spatstat]{Jest}}.
  It is hoped that this will minimise or remove any confusion for users
  of both \pkg{spatstat} and \pkg{spherstat}. 
}
\seealso{
  \code{\link[spatstat]{Jest}}, \code{\link{Ksphere}}
  \code{\link{Fsphere}}, \code{\link{Gsphere}}.
}
\references{
Baddeley, A.J. Spatial sampling and censoring. In O.E. Barndorff-Nielsen, W.S. Kendall and M.N.M. van Lieshout (eds) \emph{Stochastic Geometry: Likelihood and Computation}. Chapman and Hall, 1998. Chapter 2, pages 37--78. 

Baddeley, A.J. and Gill, R.D. The empty space hazard of a spatial pattern. Research Report 1994/3, Department of Mathematics, University of Western Australia, May 1994. 

Baddeley, A.J. and Gill, R.D. Kaplan-Meier estimators of interpoint distance distributions for spatial point processes. \emph{Annals of Statistics} \bold{25} (1997) 263--292. 

Baddeley, A., Kerscher, M., Schladitz, K. and Scott, B.T. Estimating the J function without edge correction. \emph{Statistica Neerlandica} \bold{54} (2000) 315--328. 

Borgefors, G. Distance transformations in digital images. \emph{Computer Vision, Graphics and Image Processing} \bold{34} (1986) 344--371. 

Cressie, N.A.C. \emph{Statistics for spatial data}. John Wiley and Sons, 1991. 

Diggle, P.J. \emph{Statistical analysis of spatial point patterns}. Academic Press, 1983. 

Lawrence, T.J. (2017) Master's Thesis, University of Western Australia.

Ripley, B.D. \emph{Statistical inference for spatial processes}. Cambridge University Press, 1988. 

Stoyan, D, Kendall, W.S. and Mecke, J. \emph{Stochastic geometry and its applications}. 2nd edition. Springer Verlag, 1995. 

Thonnes, E. and Van Lieshout, M.N.M, A comparative study on the power of Van Lieshout and Baddeley's J-function. \emph{Biometrical Journal} \bold{41} (1999) 721--734. 

Van Lieshout, M.N.M. and Baddeley, A.J. A nonparametric measure of spatial interaction in point patterns. \emph{Statistica Neerlandica} \bold{50} (1996) 344--361. 

}
\examples{
sph <- sphwin(type="sphere")
sph.pp <- rpoispp.sphwin(win=sph, lambda=10)
sph.ref <- rpoispp.sphwin(win=sph, lambda=150)
sph.Jest <- Jsphere(X=sph.pp, refpoints=sph.ref)
}
\keyword{spatial}
\keyword{nonparametric}

