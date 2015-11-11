\name{rHardcore.sp2}
\alias{rHardcore.sp2}
\title{
  Simulation of the Hard Core Process
}
\description{
  Generate a random pattern of points, a simulated realisation of the
  hard core process, using the Metropolis-Hastings simulation algorithm.   
}
\usage{
  rHardcore.sp2(beta, R=0, p=0.5, m=100, win = sphwin(type="sphere"), proper=TRUE)
}
\arguments{
  \item{beta}{
    intensity parameter (a positive number).
  }
  \item{R}{
    interaction parameter (a number between 0 and 1, inclusive).
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
}
\details{
  \code{rHardcore.sp2} uses the Metropolis Hasting algorithm to generate
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
\code{\link{rStrauss.sp2}}, where \code{gamma=0}, but a warning message 
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
argument \code{proper} is provided; if \code{proper=FALSE}, then
\code{m} iterations of the Metropolis-Hastings algorithm have been
completed.  If \code{proper=TRUE}, then \code{m} is ignored and as many
iterations of the Metropolis-Hastings algorithm are performed as is
needed to attain a true hard core process (i.e. with
\eqn{
  gcdist(x, X \setminus \{x\}) \ge R
}{
  gcdist(x, X \ {x}) >= R
}
for all \eqn{x \in X}{x in X}.

While it is fashionable to decry the Metropolis-Hastings algorithm for
its poor convergence and other properties, it has the advantage of being
easy to implement for a wide range of models.  
}
\value{
An object of class \code{\link{sp2}}, giving the locations of simulated points.
}
\references{
  Lawrence, T.J. (2016) Master's Thesis, University of Western Australia.
}
\author{
Tom Lawrence tjlawrence@bigpond.com
}
\note{
This function is the analogue for point processes on the sphere of the
function \code{\link{rmh}} in \code{\link{spatstat}} when used in the
third example in \code{\link{rmh.default}}, which is the corresponding
function for point processes in \code{R^2}.  Hence elements of this help
page have been taken from those for \code{\link{rmh}} and
\code{\link{rHardcore}} (also in \code{\link{spatstat}}, with the
permission of A. J. Baddeley.  This enables the code to be highly
efficient and give corresponding output to \code{\link{rmh}}, and for
the information on this help page to be consistent with that for
\code{\link{rmh}} and \code{\link{rHardcore}}.  It is hoped that this
will minimise or remove any confusion for users of both
\code{\link{spatstat}} and \code{\link{spherstat}}. 
}
\seealso{
\code{link{rHardcore}}, \code{\link{rMatClust.sp2}},
\code{\link{rMaternI.sp2}}, \code{\link{rMaternII.sp2}},
\code{\link{rpoispp.sp2}}, \code{\link{rStrauss.sp2}},
\code{\link{rThomas.sp2}} 
}
\examples{
rH1 <- rHardcore.sp2(beta=5, R=0.04*pi, p=0.5, m=100, win=sphwin())
rH1
}
\keyword{spatial}
\keyword{datagen}
