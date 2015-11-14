\name{rStrauss.sp2}
\alias{rStrauss.sp2}
\title{
  Simulation of the Strauss Process
}
\description{
  Generate a random pattern of points, a simulated realisation of the
  Strauss process, using the Metropolis-Hastings simulation algorithm.   
}
\usage{
rStrauss.sp2(beta, gamma, R=0, p=0.5, m=100, win = sphwin(type="sphere"))
}
\arguments{
  \item{beta}{
    intensity parameter (a positive number).
  }
  \item{gamma}{
    interaction strength parameter (a number between 0 and 1, inclusive).
  }
  \item{R}{
    interaction range parameter (a number between 0 and \code{pi}).
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
    f(x_1, \ldots, x_n) = alpha \, beta^{n(x)} gamma^{s(x)}
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
  pattern. If \code{gamma=1} it reduces to a Poisson process (complete
  spatial randomness) with intensity \code{beta}. If \code{gamma=0} it
  is called a ``hard core process'' with hard core radius \code{R/2},
  since no pair of points is permitted to lie closer than \code{R} units
  apart.  

  If \code{gamma=1}, a warning message will appear to advise that the
  generated process is a hard-core process, but the process well
  generate as normal.  The function \code{\link{rHardcore.sp2}}
  generates the same process but without the warning message.  It is
  recommended that \code{\link{rpoispp.sp2}} be used to generate a
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
  An object of class \code{\link{sp2}}, giving the locations of simulated points.
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
function \code{\link[spatstat]{rmh}} in \pkg{spatstat} when used in the
first example in \code{\link[spatstat]{rmh.default}}, which is the corresponding
function for point processes in \code{R^2}.  Hence elements of this help
pages have been taken from those for \code{\link[spatstat]{rmh}} and
\code{\link[spatstat]{rStrauss}} (also in \pkg{spatstat}), with the
permission of A. J. Baddeley.  This enables the code to be highly
efficient and give corresponding output to \code{\link[spatstat]{rmh}}, and for
the information on this help page to be consistent with that for
\code{\link[spatstat]{rmh}} and \code{\link[spatstat]{rStrauss}}.
It is hoped that this will minimise or remove any confusion for users of both
\pkg{spatstat} and \pkg{spherstat}. 
}
\seealso{
\code{\link{rHardcore.sp2}}, \code{\link{rMatClust.sp2}},
\code{\link{rMaternI.sp2}}, \code{\link{rMaternII.sp2}},
\code{\link{rpoispp.sp2}}, \code{\link[spatstat]{rStrauss}},
\code{\link{rThomas.sp2}} 
}
\examples{
rS1 <- rStrauss.sp2(beta=5, gamma=0.5, R=0.04*pi, p=0.5, m=100, win=sphwin())
rS1
}
\keyword{spatial}
\keyword{datagen}

