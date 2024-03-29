\name{rHardcore.sphwin}
\alias{rHardcore.sphwin}
\title{
  Simulation of the Hard Core Process
}
\description{
  Generate a random pattern of points, a simulated realisation of the
  hard core process, using the Metropolis-Hastings simulation algorithm.   
}
\usage{
  rHardcore.sphwin(beta, R=0, p=0.5, m=100, win = sphwin(type="sphere"),
 proper=TRUE, nsim=1, drop=TRUE, expand=TRUE, as.sp=TRUE, ndim="2")
}
\arguments{
  \item{beta}{
    intensity parameter (a positive number).
  }
  \item{R}{
    interaction parameter (a number between 0 and \eqn{\pi\rho}{pi*rho}, inclusive, where \eqn{\rho}{rho} is the radius of the sphere that to simulate the point pattern on.).
  }
  \item{p}{
    Metropolis-Hastings parameter: The conditional probability of
    proposing a birth given that a birth or death is proposed. 
  }
  \item{m}{
    Metropolis-Hastings parameter: The number of proposals that should
    be run. 
  }
  \item{win}{
    The window, an object of type \code{\link{sphwin}}.
  }
  \item{proper}{
    Logical.  If TRUE, then \code{m} is ignored, and the
    Metropolis-Hastings algorithm continues through as many iterations
    as is needed to achieve a hard core process.  See Details. 
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
  \code{rHardcore.sphwin} uses the Metropolis Hasting algorithm to generate
  a realisation of a hard core process inside the window \code{win}. 

The hard core process is a model for spatial inhibition, in which no
points can occur within a distance \code{R} of each other; with
interaction radius \code{R} and parameter \code{beta}, it is the
pairwise interaction point process with probability density  
\deqn{
  f(x_1, \ldots, x_n) = \alpha \, \beta^{n(x)}
}{
  f(x[1], ..., x[n]) = alpha . beta^(n(x))
}
where \eqn{x_1, \ldots, x_n}{x[1], ..., x[n]} represent the points of
the pattern, \eqn{n(x)} is the number of points in the pattern and
\eqn{\alpha}{alpha} is the normalising constant. Intuitively, each point
of the pattern contributes a factor \eqn{\beta}{beta} to the probability
density, and each pair of points closer than \eqn{r} units apart contributes a
factor \eqn{\gamma}{gamma} to the density.   

The Hardcore process may also be generated using
\code{\link{rStrauss.sphwin}}, where \code{gamma=0}, but a warning message 
is created advising that such a process is hard core rather than
Strauss. 

The Metropolis-Hastings algorithm is a Markov Chain, whose states are
spatial point patterns, and whose limiting distribution is the desired
point process. After running the algorithm for a very large number of
iterations, we may regard the state of the algorithm as a realisation
from the desired point process.  

However, there are difficulties in deciding whether the algorithm has
run for ``long enough''. The convergence of the algorithm may indeed be
extremely slow. No guarantees of convergence are given!  Hence, the
argument \code{proper} is provided; if \code{proper=FALSE}, then \code{rHardcore.sphwin} will perform 
\code{m} iterations of the Metropolis-Hastings algorithm and give the resulting point pattern.
completed.  Otherwise, (if \code{proper=TRUE}), then \code{rHardcore.sphwin} ignores \code{m} is ignored and as many
iterations of the Metropolis-Hastings algorithm are performed as is
needed to attain a true hard core process (i.e. with
\eqn{
  d(x, X \setminus \{x\}) \ge R
}{
 d(x, X \ {x}) >= R
}
for all \eqn{x \in X}{x in X}.

While it is fashionable to decry the Metropolis-Hastings algorithm for
its poor convergence and other properties, it has the advantage of being
easy to implement for a wide range of models.  
}
\value{
If \code{nsim=1} and \code{drop=FALSE} then a single item as described below; otherwise a list containing \code{nsim} items.

An item is determined by the values of \code{as.sp} and \code{ndim}:

If \code{as.sp=FALSE} and \code{ndim="2"}, a two column matrix giving the locations of the simulated points.

If \code{as.sp=FALSE} and \code{ndim="3"}, a three column matrix giving the locations of the simulated points.

If \code{as.sp=TRUE} and \code{ndim="2"}, an object of class \code{\link{sp2}} giving the locations of the simulated points.

If \code{as.sp=TRUE} and \code{ndim="3"}, an object of class \code{\link{sp3}} giving the locations of the simulated points.
}
\references{
  Lawrence, T.J. (2017) Master's Thesis, University of Western Australia.
}
\author{
Tom Lawrence tjlawrence@bigpond.com
}
\note{
This function is the analogue for point processes on the sphere of the
function \code{\link[spatstat.random]{rmh}} in \pkg{spatstat} when used in the
third example in \code{\link[spatstat.random]{rmh.default}}, which is the corresponding
function for point processes in \eqn{R^2}{R^2}.  Hence elements of this help
page have been taken from those for \code{\link[spatstat.random]{rmh}} and
\code{\link[spatstat.random]{rHardcore}} (also in \pkg{spatstat}), with the
permission of A. J. Baddeley.  This enables the code to be highly
efficient and give corresponding output to \code{\link[spatstat.random]{rmh}}, and for
the information on this help page to be consistent with that for
\code{\link[spatstat.random]{rmh}} and \code{\link[spatstat.random]{rHardcore}}.  It is hoped that this
will minimise or remove any confusion for users of both
\pkg{spatstat} and \pkg{spherstat}. 
}
\seealso{
\code{\link[spatstat.random]{rHardcore}}, \code{\link{rMatClust.sphwin}},
\code{\link{rMaternI.sphwin}}, \code{\link{rMaternII.sphwin}},
\code{\link{rpoispp.sphwin}}, \code{\link{rStrauss.sphwin}},
\code{\link{rThomas.sphwin}} 
}
\examples{
rH1 <- rHardcore.sphwin(beta=5, R=0.04*pi, p=0.5, m=100, win=sphwin())
rH1
}
\keyword{spatial}
\keyword{datagen}
