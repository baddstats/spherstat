\name{allstats.sphwin}
\alias{allstats.sphwin}
\title{
Calculate four standard summary functions of a point pattern
}
\description{
Calculates the \eqn{F}, \eqn{G}, \eqn{J}, and \eqn{K} summary functions for a point pattern on a (subset of a) sphere. Returns them as a function array (of class \code{"fasp"}, see \code{\link[spatstat]{fasp.object}}). 
}
\usage{
allstats.sphwin(X, ..., dataname=NULL, verbose=FALSE)
}
\arguments{
  \item{X}{
    The observed point pattern, from which estimates of \eqn{F(r)}, \eqn{G(r)}, \eqn{J(r)} and \eqn{K(r)} will be
    computed.  An object of class \code{"sp2"} or \code{"sp3"}, or
    a 2 or 3 column matrix giving the locations of points in spherical
    coordinates. 
  }
  \item{\dots}{
    Optional arguments, passed to the summary functions \code{\link{Fsphere}}, \code{\link{Gsphere}}, \code{\link{Jsphere}}, and \code{\link{Ksphere}}.
  }
  \item{dataname}{
    A character string giving an optional (alternative) name for the point pattern. 
  }  
  \item{verbose}{
  A logical value; if \code{TRUE}, progress reports are printed during calculation. 
  }
}
\details{
This computes four standard summary statistics for a point pattern: the empty space function \eqn{F(r)}, nearest neighbour distance distribution function \eqn{G(r)}, van Lieshout-Baddeley function \eqn{J(r)} and Ripley's function \eqn{K(r)}. The real work is done by \code{\link{Fsphere}}, \code{\link{Gsphere}}, \code{\link{Jsphere}}, and \code{\link{Ksphere}}, respectively. Consult the help files for these functions for further information about the statistical interpretation of \eqn{F}, \eqn{G}, \eqn{J} and \eqn{K}. 
If \code{verbose} is \code{TRUE}, then ``progress reports'' (just indications of completion) are printed out when the calculations are finished for each of the four function types. 
The overall title of the array of four functions (for plotting by \code{\link[spatstat]{plot.fasp}}) will be formed from the argument \code{dataname}. If this is not given, it defaults to the expression for \code{X} given in the call to \code{allstats.sphwin}. 
}
\value{
A list of length 4 containing the \eqn{F}, \eqn{G}, \eqn{J} and \eqn{K} functions respectively.
 
The list can be plotted directly using \code{plot} (which dispatches to \code{\link[spatstat]{plot.solist}}). 

Each list entry retains the format of the output of the relevant estimating routine \code{\link{Fsphere}}, \code{\link{Gsphere}}, \code{\link{Jsphere}}, or \code{\link{Ksphere}}. Thus each entry in the list is a function value table (object of class \code{"fv"}, see \code{\link[spatstat]{fv.object}}). 

}

\section{Note}{
1. Sizeable amounts of memory may be needed during the calculation

2.   This function is the analogue for point processes on the sphere of the
  function \code{\link[spatstat]{allstats}} in \code{\link{spatstat}}, which is the
  corresponding function for point processes in \eqn{R^2}{R^2}.  Indeed, almost all of the elements of 
  the code for \code{allstats.sphwin} and much of the help page have been taken from
  \code{\link[spatstat]{allstats}} with the permission of A. J. Baddeley.  This
  enables the code to be highly efficient and give corresponding output
  to, and for the information on this help page to be consistent with
  that for the function \code{\link[spatstat]{allstats}}.  It is hoped that this will
  minimise or remove any confusion for users of both
  \pkg{spatstat} and \pkg{spherstat}. 
}

\author{
  Adrian Baddeley <email: Adrian.Baddeley@curtin.edu.au>, Rolf Turner and Tom Lawrence <email:tjlawrence@bigpond.com>. 
}

\seealso{
 \code{\link{Fsphere}}, \code{\link{Gsphere}}, \code{link{Jsphere}}, \code{\link{Ksphere}}, \code{\link[spatstat]{allstats}}, \code{\link[spatstat]{plot.solist}}, \code{\link[spatstat]{plot.fv}}, \code{\link[spatstat]{fv.object}}.
}
\examples{
sph <- sphwin(type="sphere")
sph.pp <- rpoispp.sphwin(win=sph, lambda=10)
sph.all <- allstats.sphwin(sph.pp)
}
\keyword{spatial}
\keyword{nonparametric}
\keyword{hplot}
