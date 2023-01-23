\name{rStrauss.sphwin}
\alias{rStrauss.sphwin}
\title{
  Simulation of the Strauss Process
}
\description{
  Generate a random pattern of points, a simulated realisation of the
  Strauss process, using the Metropolis-Hastings simulation algorithm.   
}
\usage{
rStrauss.sphwin(beta, gamma, R=0, p=0.5, m=100, win = sphwin(type="sphere"),
 nsim=1, drop=TRUE, expand=TRUE, as.sp=TRUE, ndim="2")
}
\arguments{
  \item{beta}{
    intensity parameter (a positive number).
  }
  \item{gamma}{
    interaction strength parameter (a number between 0 and 1, inclusive).
  }
  \item{R}{
    interaction range parameter (a number between 0 and \eqn{\pi\rho}{pi*rho}, inclusive, where \eqn{\rho}{rho} is the radius of the sphere that to simulate the point pattern on).
  }
  \item{p}{
    Metropolis-Hastings parameter: The conditional probability of
    proposing a birth given that a birth or death is proposed. 
  }
  \item{m}{
    Metropolis-Hastings parameter: The number of proposals that should be run.
  }
  \item{win}{
    The window, an object of type \code{\link{sphwin}}.
  }
  \item{nsim}{
    Number of simulated realisations to be generated.
  }
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the result will be a point pattern, rather than a list containing a point pattern. 
  }
  \item{expand}{ 
    Logical. If \code{TRUE} (the default), the process is simulated on the entire sphere, and the output is only those points in the window defined by \code{win}.  If \code{FALSE}, the process is simulated within the window defined by \code{win}.
  }
  \item{as.sp}{
    Logical.  If TRUE, returns an object of class as defined by
    \code{sp.dim}.  Otherwise, returns a matrix.  See Value. 
  }
  \item{ndim}{
    A string, taking value \code{"2"} or \code{"3"}.  Specifies whether the object should contain the locations of the points in spherical coordinates (\code{ndim="2"}) or Cartesian coordinates (\code{ndim="3"}).
  }
}
\details{
  \code{rStrauss.sp2} uses the Metropolis Hastings algorithm to generate
  a realisation of a Strauss process inside the window \code{win}. 

  The Strauss process (Strauss, 1975; Kelly and Ripley, 1976) is a model
  for spatial inhibition, ranging from a strong `hard core' inhibition
  to a completely random pattern according to the value of \code{gamma}.  

  The Strauss process with interaction radius \code{R} and parameters
  \code{beta} and \code{gamma} is the pairwise interaction point process
  with probability density  
  \deqn{
    f(x_1, \ldots, x_n) = \alpha \, \beta^{n(x)} \gamma^{s(x)}
  }{
    f(x[1], ..., x[n]) = alpha * beta^(n(x)) * gamma^(s(x))
  }
  where \eqn{x_1, \ldots, x_n}{x[1], ..., x[n]} represent the points of
  the pattern, \eqn{n(x)} is the number of points in the pattern,
  \eqn{s(x)} is the number of distinct unordered pairs of points that
  are closer than \eqn{R} units apart, and \eqn{\alpha}{alpha} is the
  normalising constant. Intuitively, each point of the pattern
  contributes a factor \eqn{\beta}{beta} to the probability density, and each
  pair of points closer than \eqn{r} units apart contributes a factor
  \eqn{\gamma}{gamma} to the density.  

  The interaction parameter \code{gamma} must be less than or equal to 1
  in order that the process be well-defined (Kelly and Ripley,
  1976). This model describes an ``ordered'' or ``inhibitive''
  pattern. If \code{gamma=1} it reduces to a homogeneous spherical Poisson process (complete
  spatial randomness) with intensity \code{beta}. If \code{gamma=0} it
  is called a ``hard core process'' with hard core radius \code{R/2},
  since no pair of points is permitted to lie closer than \code{R} units
  apart.  

  If \code{gamma=1}, a warning message will appear to advise that the
  generated process is a hard-core process, but the process well
  generate as normal.  The function \code{\link{rHardcore.sphwin}}
  generates the same process but without the warning message.  It is
  recommended that \code{\link{rpoispp.sphwin}} be used to generate a
  homogeneous Poisson process in preference to this function with
  \code{gamma=0} on the basis that the former uses a more efficient
  algorithm. 

  The Metropolis-Hastings algorithm is a Markov Chain, whose states are
  spatial point patterns, and whose limiting distribution is the desired
  point process. After running the algorithm for a very large number of
  iterations, we may regard the state of the algorithm as a realisation
  from the desired point process.  

  However, there are difficulties in deciding whether the algorithm has
  run for ``long enough''. The convergence of the algorithm may indeed be
  extremely slow. No guarantees of convergence are given!  

  While it is fashionable to decry the Metropolis-Hastings algorithm for
  its poor convergence and other properties, it has the advantage of
  being easy to implement for a wide range of models.  
}
\value{
  If \code{nsim=1} and \code{drop=FALSE} then a single item as described below; otherwise a list containing \code{nsim} items.

  An item is determined by the values of \code{as.sp} and \code{ndim}:

  If \code{as.sp=FALSE} and \code{ndim="2"}, a two column matrix giving the locations of the simulated points.

  If \code{as.sp=FALSE} and \code{ndim="3"}, a three column matrix giving the locations of the simulated points.

  If \code{as.sp=TRUE} and \code{ndim="2"}, an object of class \code{sp2} giving the locations of the simulated points.

  If \code{as.sp=TRUE} and \code{ndim="3"}, an object of class \code{sp3} giving the locations of the simulated points.
}
\references{
  Kelly, F.P. and Ripley, B.D. (1976)
  On Strauss's model for clustering. \emph{Biometrika} \bold{64}, 357--360. 

  Strauss, D.J. (1975) A model for clustering. \emph{Biometrika}
  \bold{63}, 467-475.
}

\author{
Tom Lawrence <email: tjlawrence@bigpond.com>
}
\note{
This function is the analogue for point processes on the sphere of the
function \code{\link[spatstat.random]{rmh}} in \pkg{spatstat} when used in the
first example in \code{\link[spatstat.random]{rmh.default}}, which is the corresponding
function for point processes in \eqn{R^2}{R^2}.  Hence elements of this help
pages have been taken from those for \code{\link[spatstat.random]{rmh}} and
\code{\link[spatstat.random]{rStrauss}} (also in \pkg{spatstat}), with the
permission of A. J. Baddeley.  This enables the code to be highly
efficient and give corresponding output to \code{\link[spatstat.random]{rmh}}, and for
the information on this help page to be consistent with that for
\code{\link[spatstat.random]{rmh}} and \code{\link[spatstat.random]{rStrauss}}.
It is hoped that this will minimise or remove any confusion for users of both
\pkg{spatstat} and \pkg{spherstat}. 
}
\seealso{
\code{\link{rHardcore.sphwin}}, \code{\link{rMatClust.sphwin}},
\code{\link{rMaternI.sphwin}}, \code{\link{rMaternII.sphwin}},
\code{\link{rpoispp.sphwin}}, \code{\link[spatstat.random]{rStrauss}},
\code{\link{rThomas.sphwin}} 
}
\examples{
rS1 <- rStrauss.sphwin(beta=5, gamma=0.5, R=0.04*pi, p=0.5, m=100, win=sphwin())
rS1
}
\keyword{spatial}
\keyword{datagen}

