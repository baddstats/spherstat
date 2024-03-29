\name{Ksphere}
\alias{Ksphere}
\title{
  Estimate of Spherical K Function
}
\description{
  Given a point pattern on the sphere, or on a region of the sphere,
  compute an estimate of the \eqn{K}-function.
}
\usage{
  Ksphere(X, win=sphwin(), r,
          correction=c("un", "iso", "rs", "rsm"),
          ratio=FALSE, lambdavalues=NULL, update=TRUE)
}
\arguments{
  \item{X}{
    The observed point pattern, from which an estimate of \eqn{K(r)} will be
    computed.  An object of class \code{sp2} or \code{sp3}, or
    a 2 or 3 column matrix giving the locations of points in spherical
    coordinates. 
  }
  \item{win}{
    The window; an object of type \code{sphwin}.  Only required if
    \code{X} is an object of class \code{matrix}. 
  }
  \item{r}{
    Optional.  Numeric vector. The values of the argument \eqn{r} at
    which \eqn{K(r)} should be evaluated.  
  }
  \item{correction}{
    Character string or character vector specifying the choice of
    edge corrections.
  }
  \item{ratio}{
    Logical value specifying whether to store information about the
    numerator and denominator of the estimator, for use by
    other functions.
  }
  \item{lambdavalues}{
    Optional.  If this argument is given, the inhomogeneous
    \eqn{K} function is computed.
    Either a numeric vector containing estimates of the intensity
    at each data point, or a fitted point process model (object of class
    \code{"sphppm"}) which predicts the intensity.
  }
  \item{update}{
    Logical value indicating, when \code{lambdavalues} is a fitted
    point process model, whether to re-fit the model to
    the point pattern \code{X} before extracting the fitted intensity
    values. 
  }
}
\details{
If \code{lambdavalues=NULL}, this computes the homogeneous \eqn{K} function for a point pattern on the sphere, analogous to that for \eqn{R^2}{R^2} proposed by Ripley (1976).

If \code{lambdavalues!=NULL}, this computes a generalisation of the \eqn{K} function for inhomogeneous point patterns, analagous to that proposed by Baddeley, Moller and Waagepetersen (2000). 

The homogeneous or Ripley's \eqn{K} function (variously known as the "reduced second-order moment function" and "Ripley's \eqn{K} function" of a stationary point process \code{X} is defined so that \eqn{\lambda K(r)}{lambda K(r)} equals the expected number of additional random points within a distance \eqn{r} of a typical random point of \code{X}. Here \eqn{\lambda}{lambda} is the intensity of the process, i.e. the expected number of points of \code{X} per unit area. The \eqn{K} function is determined by the second order moment properties of \code{X}. 

An estimate of \eqn{K} derived from a spatial point pattern dataset can be used in exploratory data analysis and formal inference about the pattern (Cressie, 1991; Diggle, 1983; Ripley, 1977, 1988). In exploratory analyses, the estimate of \eqn{K} is a useful statistic summarising aspects of inter-point "dependence" and "clustering". For inferential purposes, the estimate of \eqn{K} is usually compared to the true value of \eqn{K} for a completely random (Poisson) point process, which is 

\deqn{K(r) = 2\pi\rho^{2}(1 - \cos(r/\rho))}{K(r) = 2 * pi * rho^2 * (1 - cos(r/rho))}

where \eqn{\rho}{rho} is the radius of the sphere on which \code{X} is observed. Deviations between the empirical and theoretical \eqn{K} curves may suggest spatial clustering or spatial regularity. 

Provided that \code{lambdavalues=NULL}, this routine \code{Ksphere} estimates the \eqn{K} function of a homogeneous point process, given observation of the process inside a window that is a sphere, cap, band, band complement, wedge, or polygon where all boundaries are great circle arcs. The argument \code{X} is interpreted as a point pattern object (of class \code{"sp2"} or \code{"sp3"}; or of class \code{"matrix"} with 2 or 3 columns giving the coordinates of points, and accompanied by an object of class \code{"sphwin"} giving the observation window). 

The estimation of \eqn{K} is hampered by edge effects arising from the unobservability of points of the random pattern outside the window. An edge correction is needed to reduce bias (Baddeley, 1998; Ripley, 1988). The corrections implemented here are 

\code{rs}, \code{rsm}
the border and modified border methods or "reduced sample" estimators (see Ripley, 1988). These is the least efficient (statistically) and the fastest to compute.

\code{iso}
The analogue of Ripley's isotropic correction (see Ripley, 1988; Ohser, 1983).  This is implemented for caps, bands, band complements, wedges and polygons where all boundaries are great circle arcs.

\code{un}
Uncorrected estimate. An estimate of the \eqn{K} function without edge correction. (i.e. setting \eqn{e[i,j] = 1}{e[i,j] = 1} in the equation below. Unless the window is a sphere, this estimate is biased and should not be used for data analysis, unless you have an extremely large point pattern (more than 100,000 points). 

The estimates of \eqn{K(r)} are of the form 
\deqn{
  Khom(r) = (a/(n * (n-1))) * \sum[i,j] I(d[i,j] <= r) e[i,j])
}{
  Khom(r) = (a/(n * (n-1))) * sum[i,j] I(d[i,j] <= r) e[i,j])
}

where \eqn{a}{a} is the area of the window, \eqn{n}{n} is the number of data points, and the sum is taken over all ordered pairs of points \eqn{i}{i} and \eqn{j}{j} in \code{X}. Here \eqn{d[i,j]}{d[i,j]} is the distance between the two points, and \eqn{I(d[i,j] <= r)}{I(d[i,j] <= r)} is the indicator that equals \code{1} if the distance is less than or equal to \eqn{r}{r}. The term \eqn{e[i,j]}{e[i,j]} is the edge correction weight (which depends on the choice of edge correction listed above). 


The inhomogeneous \eqn{K} function \eqn{Kinhom(r)} is a direct generalisation to nonstationary point processes. Suppose \code{x} is a point process with non-constant intensity \eqn{\lambda(u)}{lambda(u)} at each location \code{u}. Define \eqn{Kinhom(r)} to be the expected value, given that \code{u} is a point of \code{x}, of the sum of all terms \eqn{1/\lambda(x[j])}{1/lambda(x[j])} over all points \code{x[j]} in the process separated from \code{u} by a distance less than \code{r}. This reduces to the ordinary \eqn{K} function if \eqn{\lambda(u)}{lambda()} is constant. If \code{x} is an inhomogeneous Poisson process with intensity function \eqn{\lambda(u)}{lambda(u)}, then \eqn{Kinhom(r) = 2\pi\rho^{2}(1-\cos(r))}{Kinhom(r) = 2 * pi * rho^2 * (1-cos(r))}

Given a point pattern dataset, the inhomogeneous \eqn{K} function can be estimated essentially by summing the values \eqn{1/(\lambda(x[i]) * \lambda(x[j]))}{1/(lambda(x[i]) * lambda(x[j]))} for all pairs of points \eqn{(x[i],x[j])}{(x[i],x[j])} separated by a distance less than \eqn{r}{r}. 

This allows us to inspect a point pattern for evidence of interpoint interactions after allowing for spatial inhomogeneity of the pattern. Values \eqn{Kinhom(r) > 2\pi\rho^{2} (1-\cos(r/\rho))}{Kinhom(r) > 2 * pi * rho^2 * (1-cos(r/\rho))} are suggestive of clustering. 

The argument \code{lambdavalues} should supply the (estimated) values of the intensity function \code{lambda}. It may be either 
a numeric vector containing the values of the intensity function at the points of the pattern \code{X}. 
a fitted point process model
(object of class "\code{\link{sphppm}}") whose fitted trend can be used as the fitted intensity.

If \code{lambdavalues} is a numeric vector, then its length should be equal to the number of points in the pattern \code{X}. The value \code{lambda[i]} is assumed to be the the (estimated) value of the intensity \eqn{\lambda (x_i)}{lambda(x_i)} for the point \code{x[i]} of the pattern \code{X}. Each value must be a positive number; \code{NA}'s are not allowed. 

Edge corrections are used to correct bias in the estimation of \eqn{K}. Each edge-corrected estimate of \eqn{K(r)} is of the form 
\deqn{
  Kinhom(r) = \sum_{i} \sum_{j} 1(d[i,j] <= r) * e(x[i],x[j],r)/(\lambda(x[i]) * \lambda(x[j]))
}{
  Kinhom(r) = sum_i sum_j 1(d[i,j] <= r) * e(x[i],x[j],r)/(lambda(x[i]) * lambda(x[j]))
}

where \eqn{d[i,j]}{d[i,j]} is the minimum great circle distance between points \eqn{x[i]}{x[i]} and \eqn{x[j]}{x[j]}, and \eqn{e(x[i],x[j],r)}{e(x[i],x[j],r)} is an edge correction factor.  The only available correction is the 'isotropic' correction, for which

\deqn{
  e(x[i],x[j],r) =  1/(area(W) g(x[i],x[j]))
}{
  e(x[i],x[j],r) =  1/(area(W) g(x[i],x[j]))
}

where \eqn{g(x[i],x[j])}{g(x[i],x[j])} is the fraction of the circumference of the circle with centre \eqn{x[i]}{x[i]} and radius \eqn{d(x[i],x[j])}{d(x[i],x[j])} which lies inside the window. 

}
\value{
  An object of class \code{"fv"}, see \code{\link[spatstat.explore]{fv.object}},
  which can be plotted directly using \code{\link[spatstat.explore]{plot.fv}}. 
}
\author{
  Tom Lawrence.
}
\seealso{
  \code{\link{Fsphere}}, \code{\link{Gsphere}}, \code{\link{Jsphere}}, \code{\link[spatstat.explore]{Kest}}, \code{\link[spatstat.explore]{Kinhom}}
}
\references{
Baddeley, A.J. Spatial sampling and censoring. In O.E. Barndorff-Nielsen, W.S. Kendall and M.N.M. van Lieshout (eds) \emph{Stochastic Geometry: Likelihood and Computation}. Chapman and Hall, 1998. Chapter 2, pages 37--78. 

Baddeley, A., Moller, J. and Waagepetersen, R. (2000) Non- and semiparametric estimation of interaction in inhomogeneous point patterns. \emph{Statistica Neerlandica} \bold{54}, 329--350. 

Cressie, N.A.C. \emph{Statistics for spatial data}. John Wiley and Sons, 1991. 

Diggle, P.J. \emph{Statistical analysis of spatial point patterns}. Academic Press, 1983. 

Lawrence, T. (2017) Master's Thesis, University of Western Australia.

Ohser, J. (1983) On estimators for the reduced second moment measure of point processes. \emph{Mathematische Operationsforschung und Statistik, series Statistics}, \bold{14}, 63 -- 71. 

Ohser, J. and Stoyan, D. (1981) On the second-order and orientation analysis of planar stationary point processes. \emph{Biometrical Journal} \bold{23}, 523--533. 

Raskin, R.G. (1994) Spatial analysis on the sphere: a review (94-7) \emph{NGCIA technical report} (Available from http://eprints.cdlib.org/uc/item/5748n2xz).  Accessed 4 November 2015.

Ripley, B.D. (1977) Modelling spatial patterns (with discussion). \emph{Journal of the Royal Statistical Society}, Series B, \bold{39}, 172 -- 212. 

Ripley, B.D. \emph{Statistical inference for spatial processes}. Cambridge University Press, 1988. 

Stoyan, D, Kendall, W.S. and Mecke, J. (1995) \emph{Stochastic geometry and its applications}. 2nd edition. Springer Verlag. 

Stoyan, D. and Stoyan, H. (1994) \emph{Fractals, random shapes and point fields: methods of geometrical statistics}. John Wiley and Sons. 
}
\note{
  This function is the analogue for point processes on the sphere of the
  functions \code{\link[spatstat.explore]{Kest}} and \code{\link[spatstat.explore]{Kinhom}} in \pkg{spatstat}, which is the
  corresponding function for point processes in \eqn{R^2}{R^2}.  Hence elements of
  the code for \code{Ksphere} and help page have been taken from
  \code{\link[spatstat.explore]{Kest}} and \code{\link[spatstat.explore]{Kinhom}} with the permission of A. J. Baddeley.  This
  enables the code to be highly efficient and give corresponding output
  to, and for the information on this help page to be consistent with
  that for the functions \code{\link[spatstat.explore]{Kest}} and \code{\link[spatstat.explore]{Kinhom}}.  It is hoped that this will
  minimise or remove any confusion for users of both
  \pkg{spatstat} and \pkg{spherstat}. 
}
\examples{
   win <- sphwin()
   X <- runif.sphwin(100, win)
   rv <- seq(0, pi/2, length=512)
   K <- Ksphere(X, win, rv)
   K
   plot(K)

   sphcap <- sphwin(type="band", param=c(0, pi/3), ref=c(0,0))
   Xcap <- runif.sphwin(100, sphcap)
   rv3 <- seq(0, pi/3, length=512)
   Kcap <- Ksphere(Xcap, r=rv3)
   plot(Kcap, xlim=c(0, pi/4))

   sphband <- sphwin(type="band", param=c(pi/3, 2*pi/3), ref=c(0,0))
   Xb <- runif.sphwin(100, sphband)
   Kb <- Ksphere(Xb, r=rv3)
   plot(Kb)

   sphbc <- sphwin(type="bandcomp", param=c(pi/3, 2*pi/3), ref=c(0,0))
   Xbc <- runif.sphwin(100, sphbc)
   Kbc <- Ksphere(Xbc, r=rv3)
   plot(Kbc)
}
\keyword{spatial}
\keyword{nonparametric}
