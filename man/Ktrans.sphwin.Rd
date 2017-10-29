\name{Ktrans}
\alias{Kstab.sphwin}
\alias{Kinv.sphwin}
\title{
 Transformations of Spherical K Function
}
\description{
  Given a point pattern on the sphere, or on a region of the sphere,
  compute a transform of an estimate of the \eqn{K}-function.
}
\usage{
  Kstab.sphwin(X, ...)
  Kinv.sphwin(X, ...)
}
\arguments{
  \item{X}{
    The observed point pattern, from which an estimate of \eqn{K(r)} will be
    computed.  An object of class \code{\link{sp2}} or \code{\link{sp3}}, or
    a 2 or 3 column matrix giving the locations of points in spherical
    coordinates. 
  }
  \item{...}{
    Other arguments passed to \code{\link{Ksphere}}.
  }
}
\details{
This command computes an estimate of the variance-stabilising and inverse transformations of Ripley's (1977) \eqn{K} function for the spatial point pattern X. The variance-stabilising transformation of Ripley's K-function is
\deqn{
Kstab(r) = \sqrt{(K(r)/(2*\pi*\rho^2))}
}{
Kstab(r) = sqrt(K(r)/(2*pi*rho^2))
}
while the inverse transformation is
\deqn{
Kinv(r) = \rho*\arccos(K(r)/(2*\pi*\rho^2))
}{
Kinv(r) = rho*arccos(K(r)/(2*pi*rho^2))
}
where \eqn{K(r)} is the \eqn{K}-function and \eqn{\rho}{rho} is the radius of the sphere on which \code{X} is observed. 

See \code{\link{Ksphere}} for information about Ripley's \eqn{K}-function. The variance-stabilising transformation was proposed by Besag (1977); the inverse transformation by Lawrence (2016). 

The commands \code{Kstab.sphwin} and \code{Kinv.sphwin} first calls \code{\link{Ksphere}} to compute the estimate of the \eqn{K}-function, and then applies the appropriate transformation. 

For a completely random (uniform Poisson) point pattern, the theoretical value of the inverse transformation function is \eqn{Kinv(r) = r}. The variance-stabilising transformation is more appropriate for use in simulation envelopes and hypothesis tests. 


Other plots that can be useful include centring (include the argument \code{.-theo~r}) the function, or for either of \eqn{K(r)} or \eqn{Kstab(r)}, plotting the estimated function against the theoretical value for a model (use e.g. the argument \code{.~theo} in the plot function).

See \code{\link{Ksphere}} for the list of arguments. 
}

\value{
  An object of class \code{"fv"}, see \code{\link[spatstat]{fv.object}},
  which can be plotted directly using \code{\link[spatstat]{plot.fv}}. 
}
\author{
  Tom Lawrence.
}
\seealso{
  \code{\link{Fsphere}}, \code{\link{Gsphere}}, \code{\link{Jsphere}}, \code{\link[spatstat]{Kest}}, \code{\link[spatstat]{Kinhom}}, \code{\link{Fstab.sphwin}}, \code{\link{Finv.sphwin}}, \code{\link{Gstab.sphwin}}, \code{\link{Ginv.sphwin}}, \code{\link[spatstat]{Lest}}
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
  functions \code{\link[spatstat]{Kest}} and \code{\link[spatstat]{Kinhom}} in \pkg{spatstat}, which is the
  corresponding function for point processes in \eqn{R^2}{R^2}.  Hence elements of
  the code for \code{Ksphere} and help page have been taken from
  \code{\link[spatstat]{Kest}} and \code{\link[spatstat]{Kinhom}} with the permission of A. J. Baddeley.  This
  enables the code to be highly efficient and give corresponding output
  to, and for the information on this help page to be consistent with
  that for the functions \code{\link[spatstat]{Kest}} and \code{\link[spatstat]{Kinhom}}.  It is hoped that this will
  minimise or remove any confusion for users of both
  \pkg{spatstat} and \pkg{spherstat}. 
}
\examples{
   win <- sphwin()
   X <- runif.sphwin(100, win)
   rv <- seq(0, pi/2, length=512)
   K <- Ksphere(X, win, rv)
   Kstab <- Kstab.sphwin(X, win, rv)
   Kinv <- Kinv.sphwin(X, win, rv)
   plot(K)
   plot(K, .~theo)
   plot(K, .-theo~r)
   plot(Kstab)
   plot(Kstab, .~theo)
   plot(Kinv) 
}
\keyword{spatial}
\keyword{nonparametric}
