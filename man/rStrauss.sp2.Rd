\name{rStrauss.sp2}
\alias{rStrauss.sp2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of the Strauss Process
}
\description{
Generate a random pattern of points, a simulated realisation of the Strauss process, using the Metropolis-Hastings simulation algorithm.  
}
\usage{
rStrauss.sp2(beta, gamma, R=0, p=0.5, m=100, win = sphwin())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
intensity parameter (a positive number).
}
  \item{gamma}{
interaction parameter (a number between 0 and 1, inclusive).
}
  \item{R}{
interaction parameter (a number between 0 and 1, inclusive).
}
  \item{p}{
Metropolis-Hastings parameter: The conditional probability of proposing a birth given that a birth or death is proposed.
}
  \item{m}{
Metropolis-Hastings parameter: The number of proposals that should be run.
}
  \item{win}{
The window, an object of type \code{\link{sphwin}}.
}
}
\details{
\code{rStrauss.sp2} uses the Metropolis Hasting algorithm to generate a realisation of a Strauss process inside the window \code{win}.

The Strauss process (Strauss, 1975; Kelly and Ripley, 1976) is a model for spatial inhibition, ranging from a strong `hard core' inhibition to a completely random pattern according to the value of \code{gamma}. 

The Strauss process with interaction radius \code{R} and parameters \code{beta} and \code{gamma} is the pairwise interaction point process with probability density 

\code{ f(x_1, ..., x_n) = alpha . beta^n(x) gamma^s(x) }

where \code{x_1, ..., x_n} represent the points of the pattern, \code{n(x)} is the number of points in the pattern, \code{s(x)} is the number of distinct unordered pairs of points that are closer than \code{R} units apart, and \code{alpha} is the normalising constant. Intuitively, each point of the pattern contributes a factor \code{beta} to the probability density, and each pair of points closer than \code{r} units apart contributes a factor \code{gamma} to the density. 

The interaction parameter gamma must be less than or equal to 1 in order that the process be well-defined (Kelly and Ripley, 1976). This model describes an ``ordered'' or ``inhibitive'' pattern. If \code{gamma=1} it reduces to a Poisson process (complete spatial randomness) with intensity \code{beta}. If \code{gamma=0} it is called a ``hard core process'' with hard core radius \code{R/2}, since no pair of points is permitted to lie closer than \code{R} units apart. 

If \code{gamma=1}, a warning message will appear to advise that the generated process is a hard-core process, but the process well generate as normal.  The function \code{\link{rHardcore.sp2}} generates the same process but without the warning message.  It is recommended that \code{\link{rpoispp.sp2}} be used to generate a homogeneous Poisson process in preference to this function with \code{gamma=0} on the basis that the former uses a more efficient algorithm.

The Metropolis-Hastings algorithm is a Markov Chain, whose states are spatial point patterns, and whose limiting distribution is the desired point process. After running the algorithm for a very large number of iterations, we may regard the state of the algorithm as a realisation from the desired point process. 

However, there are difficulties in deciding whether the algorithm has run for ``long enough''. The convergence of the algorithm may indeed be extremely slow. No guarantees of convergence are given! 

While it is fashionable to decry the Metropolis-Hastings algorithm for its poor convergence and other properties, it has the advantage of being easy to implement for a wide range of models. 

}
\value{
An object of class \code{\link{sp2}}, giving the locations of simulated points.
}
\references{
Kelly, F.P. and Ripley, B.D. (1976) On Strauss's model for clustering. Biometrika 63, 357-360. 
Strauss, D.J. (1975) A model for clustering. Biometrika 63, 467-475.
}
\author{
Tom Lawrence tjlawrence@bigpond.com
}
\note{
This function is the analogue for point processes on the sphere of the function \code{\link{rmh}} in \code{\link{spatstat}} when used in the first example in \code{\link{rmh.default}}, which is the corresponding function for point processes in \code{R^2}.  Hence elements of this help pages have been taken from those for \code{\link{rmh}} and \code{\link{rStrauss}} (also in \code{\link{spatstat}}, with the permission of A. J. Baddeley.  This enables the code to be highly efficient and give corresponding output to \code{\link{rmh}}, and for the information on this help page to be consistent with that for \code{\link{rmh}} and \code{\link{rStrauss}}.  It is hoped that this will minimise or remove any confusion for users of both \code{\link{spatstat}} and \code{\link{spherstat}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rHardcore.sp2}}, \code{\link{rMatClust.sp2}}, \code{\link{rMaternI.sp2}}, \code{\link{rMaternII.sp2}}, \code{\link{rpoispp.sp2}}, \code{\link{rStrauss}}, \code{\link{rThomas.sp2}}
}
\examples{
rS1 <- rStrauss.sp2(beta=50, gamma=5, R=0.04*pi, p=0.5, m=100, win=sphwin())
rS1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
