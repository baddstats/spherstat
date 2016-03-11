\name{Gtrans}
\alias{Ginv.sphwin}
\alias{Gstab.sphwin}
\title{
Transforms of the Empty Space Function F 
}
\description{
Estimates the summary function \code{G(r)} for a point pattern in a window of arbitrary shape on a (subset of) a sphere and performs either the variance-stabilising or inverse transform.
}
\usage{
Ginv.sphwin(X, ...)

Gstab.sphwin(X, ...)
}
\arguments{
  \item{X}{
    The observed point pattern, from which an estimate of \eqn{G(r)} will be
    computed and then transformed.  An object of class \code{"sp2"} or \code{"sp3"}, or
    a 2 or 3 column matrix giving the locations of points in spherical
    coordinates. 
  }
  \item{...}{
    Other arguments required by the function \code{\link{Gsphere}}.
  }
}
\details{
  The nearest neighbour distance distribution function (also called the "event-to-event", or the "inter-event" distribution) of a point process \code{X} is the cumulative distribution function \code{G} of the distance from a typical random point of \code{X} to the nearest other point of \code{X}. 

  The variance-stabilising transform of \code{G} is \code{asin(sqrt(G(r)))}.

  Alternatively, \code{Ginv.sphwin} calculates the inverse of the estimated \eqn{G} function, under the assumption that the underlying process is Poisson.  The calculated function is therefore independent of the estimated intensity of the point pattern and therefore any deviation test performed using this estimate will avoid the related conservatism issue.  The plot is also a P-P plot.

  Plots that can be useful include centring (include the argument \code{.-theo~r}) the function, or for either of \eqn{G(r)} or \eqn{Gstab(r)}, plotting the estimated function against the theoretical value for a model (use e.g. the argument \code{.~theo} in the plot function).
}
\value{
  An object of class \code{"fv"}, see \code{\link[spatstat]{fv.object}},
  which can be plotted directly using \code{\link[spatstat]{plot.fv}}. 

  Essentially a data frame containing some or all of the following columns: 
  \item{r}{the values of the argument r at which the transformation of\eqn{G(r)}
  has been estimated} 
 
  \item{rs}{the transformation of the "reduced sample" (border correction) estimator of \eqn{G(r)}}
 
  \item{km}{the trasnformation of the spatial Kaplan-Meier estimator of \eqn{G(r)}}
 
  \item{hazard}{the transformation of the hazard rate \eqn{lambda(r)} of \eqn{G(r)} by the spatial Kaplan-Meier method}
 
  \item{raw}{the transform of the uncorrected estimate of \eqn{G(r)}, , i.e. the empirical distribution of the distances from each point in the pattern X to the nearest other point of the pattern}
 
\item{han}{the transform of the Hanisch correction estimator of \eqn{G(r)}}
 
\item{theo}{the transform of the theoretical value of \eqn{G(r)} for a stationary Poisson process of the same estimated intensity.}
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

\seealso{
  \code{\link{Gsphere}}, \code{link{Jsphere}}, \code{\link{Ksphere}}, \code{\link[spatstat]{Gest}}, \code{\link{Fstab.sphwin}}, \code{\link{Finv.sphwin}}, \code{\link{Kstab.sphwin}}, \code{\link{Kinv.sphwin}}
}
\examples{
sph <- sphwin(type="sphere")
sph.pp <- rpoispp.sphwin(win=sph, lambda=10)
sph.ref <- rpoispp.sphwin(win=sph, lambda=150)
sph.Gstab <- Gstab.sphwin(X=sph.pp)
sph.Ginv <- Ginv.sphwin(X=sph.pp)
}
\keyword{spatial}
\keyword{nonparametric}
