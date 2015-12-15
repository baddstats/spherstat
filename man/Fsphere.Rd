\name{Fsphere}
\alias{Fsphere}
\alias{Fstab.sph}
\title{
The Empty Space Function F 
}
\description{
Estimates the summary function \code{F(r)} for a point pattern in a window of arbitrary shape on a (subset of) a sphere.
}
\usage{
Fsphere(X, refpoints, win = sphwin(type = "sphere"), r = NULL, ...)

Fstab.sph(X, refpoints, win = sphwin(type = "sphere"), r = NULL, ...)
}
\arguments{
  \item{X}{
    The observed point pattern, from which an estimate of \eqn{F(r)} will be
    computed.  An object of class \code{"sp2"} or \code{"sp3"}, or
    a 2 or 3 column matrix giving the locations of points in spherical
    coordinates. 
  }
  \item{refpoints}{
    The reference points in the window used to estimate the \eqn{F}
    function.  An object of class \code{"sp2"} or \code{"sp3"},
    or a 2 or 3 column matrix giving the locations of
    points in spherical coordinates. 
  }  
  \item{win}{
    The window; an object of class \code{"sphwin"}.  Only required if
    \code{X} is an object of class \code{matrix}. 
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
}
\details{
  The empty space function (also called the
  "spherical contact distribution" or the "point-to-nearest-event"
  distribution) of a stationary point process \code{X} is the cumulative
  distribution function \eqn{F} of the distance from a fixed point in
  space to the nearest point of \code{X}.  

  An estimate of \eqn{F} derived from a spatial point pattern dataset can
  be used in exploratory data analysis and formal inference about the
  pattern (Cressie, 1991; Diggle, 1983; Ripley, 1988). In exploratory
  analyses, the estimate of \code{F} is a useful statistic summarising the
  sizes of gaps in the pattern. For inferential purposes, the estimate of
  \code{F} is usually compared to the true value of \code{F} for a
  completely random (Poisson) point process, which is  
  \deqn{
    F(r)=1-\exp(-2\lambda\pi(1-\cos(r/rad)))
  }{
    F(r) = 1 - exp( - 2 * lambda * pi * (1 - cos (r/rad) ) )
  }
  where \eqn{\lambda}{lambda} is the intensity (expected number of
  points per unit area). Deviations between the empirical and
  theoretical \eqn{F} curves may suggest spatial clustering or spatial
  regularity.  

  This algorithm estimates the empty space function \eqn{F} from the
  point pattern \code{X}. It assumes that \code{X} can be treated as a
  realisation of a stationary (spatially homogeneous) random spatial
  point process in the plane, observed through a bounded window. The
  window (which is specified in \code{X}) may have arbitrary shape.  

  The argument \code{X} is interpreted as a point pattern object (of
  class \code{"ppp"}, see \code{\link[spatstat]{ppp.object}}) and can be
  supplied in any of the formats recognised by \code{\link[spatstat]{as.ppp}}.  

  The algorithm uses a discrete approximation which is controlled by
  the spacing of values of \code{r}. (See below for details.) First-time
  users are strongly advised not to specify these arguments.  

  The estimation of \eqn{F} is hampered by edge effects arising from the
  unobservability of points of the random pattern outside the window. An
  edge correction is needed to reduce bias (Baddeley, 1998; Ripley,
  1988). The edge corrections implemented here are the border method or
  "reduced sample" estimator, the spatial Kaplan-Meier estimator
  (Baddeley and Gill, 1997) and the Chiu-Stoyan estimator (Chiu and
  Stoyan, 1998).  

  The argument \code{r} is the vector of values for the distance \eqn{r}
  at which \eqn{F(r)} should be evaluated. It is also used to determine
  the breakpoints (in the sense of \code{hist}) for the computation of
  histograms of distances. The estimators are computed from histogram
  counts. This introduces a discretisation error which is controlled by
  the fineness of the breakpoints.  

  First-time users would be strongly advised not to specify
  \code{r}. However, if it is specified, \code{r} must satisfy
  \code{r[1] = 0}, and \code{max(r)} must be larger than the radius of
  the largest disc contained in the window. Furthermore, the spacing of
  successive \code{r} values must be very fine.  

  The algorithm also returns an estimate of the hazard rate function,
  \code{lambda(r)}, of \code{F(r)}. The hazard rate is defined by  
  \deqn{
    \lambda(r) = - \frac{d}{dr} \log(1 - F(r))
  }{
    lambda(r) = - (d/dr) log(1 - F(r))
  }

  The hazard rate of \eqn{F} has been proposed as a useful exploratory
  statistic (Baddeley and Gill, 1994). The estimate of \eqn{lambda(r)} given
  here is a discrete approximation to the hazard rate of the
  Kaplan-Meier estimator of \eqn{F}. Note that \eqn{F} is absolutely
  continuous (for any stationary point process \code{X}), so the hazard
  function always exists (Baddeley and Gill, 1997).  

  The naive empirical distribution of distances from each location in the
  window to the nearest point of the data pattern, is a biased estimate of
  \eqn{F}. However this is also returned by the algorithm (if
  \code{correction="none"}), as it is sometimes useful in other
  contexts. Care should be taken not to use the uncorrected empirical
  \eqn{F} as if it were an unbiased estimator of \eqn{F}.  
}
\value{
  An object of class \code{"fv"}, see \code{\link[spatstat]{fv.object}},
  which can be plotted directly using \code{\link[spatstat]{plot.fv}}. 

  Essentially a data frame containing some or all of the following columns: 
  \item{r}{the values of the argument r at which the function \eqn{F(r)}
  has been estimated} 
 
  \item{rs}{the "reduced sample" (border correction) estimator of \eqn{F(r)}}
 
  \item{km}{the spatial Kaplan-Meier estimator of \eqn{F(r)}}
 
  \item{hazard}{the hazard rate \eqn{lambda(r)} of \eqn{F(r)} by the spatial Kaplan-Meier method}
 
  \item{raw}{the uncorrected estimate of \eqn{F(r)}, i.e. the empirical distribution of the distance from a fixed point in the window to the nearest point of \code{X}}
 
\item{han}{the Hanisch correction estimator of \eqn{F(r)}}
 
\item{theo}{the theoretical value of \eqn{F(r)} for a stationary Poisson process of the same estimated intensity.}
}

\section{Note}{
Sizeable amounts of memory may be needed during the calculation
}

\references{
Baddeley, A.J. Spatial sampling and censoring. In O.E. Barndorff-Nielsen, W.S. Kendall and M.N.M. van Lieshout (eds) \emph{Stochastic Geometry: Likelihood and Computation}. Chapman and Hall, 1998. Chapter 2, pages 37--78. 

Baddeley, A.J. and Gill, R.D. The empty space hazard of a spatial pattern. Research Report 1994/3, Department of Mathematics, University of Western Australia, May 1994. 

Baddeley, A.J. and Gill, R.D. Kaplan-Meier estimators of interpoint distance distributions for spatial point processes. \emph{Annals of Statistics} \bold{25} (1997) 263--292. 

Borgefors, G. Distance transformations in digital images. \emph{Computer Vision, Graphics and Image Processing} \bold{34} (1986) 344--371. 

Chiu, S.N. and Stoyan, D. (1998) Estimators of distance distributions for spatial patterns. \emph{Statistica Neerlandica} \bold{52}, 239--246. 

Cressie, N.A.C. \emph{Statistics for spatial data}. John Wiley and Sons, 1991. 

Diggle, P.J. \emph{Statistical analysis of spatial point patterns}. Academic Press, 1983. 

Lawrence, T.J. (2016) Master's Thesis, University of Western Australia.

Ripley, B.D. \emph{Statistical inference for spatial processes}. Cambridge University Press, 1988. 

Stoyan, D, Kendall, W.S. and Mecke, J. \emph{Stochastic geometry and its applications}. 2nd edition. Springer Verlag, 1995. 

}
\author{
  Tom Lawrence <email:tjlawrence@bigpond.com>
  and Adrian Baddeley <email: Adrian.Baddeley@curtin.edu.au>
}
\note{
  This function is the analogue for point processes on the sphere of the
  function \code{\link[spatstat]{Fest}} in \code{\link{spatstat}}, which is the
  corresponding function for point processes in \eqn{R^2}{R^2}.  Hence elements of
  the code for \code{Fsphere} and help page have been taken from
  \code{\link[spatstat]{Fest}} with the permission of A. J. Baddeley.  This
  enables the code to be highly efficient and give corresponding output
  to, and for the information on this help page to be consistent with
  that for the function \code{\link[spatstat]{Fest}}.  It is hoped that this will
  minimise or remove any confusion for users of both
  \pkg{spatstat} and \pkg{spherstat}. 
}
\seealso{
  \code{\link{Gsphere}}, \code{link{Jsphere}}, \code{\link{Ksphere}}, \code{\link[spatstat]{Gest}}
}
\examples{
sph <- sphwin(type="sphere")
sph.pp <- rpoispp.sphwin(win=sph, lambda=10)
sph.ref <- rpoispp.sphwin(win=sph, lambda=150)
sph.Fest <- Fsphere(X=sph.pp, refpoints=sph.ref)
}
\keyword{spatial}
\keyword{nonparametric}
