\name{splitcaps}
\alias{splitcaps}
\title{
  Extract Points Falling in a Given Cap
}
\description{
  Extract all the points in a band (two disjoint parallel caps) that fall
  in a specified cap. 
}
\usage{
splitcaps(X, win=sphwin(type="sphere"), upper=TRUE, lower=TRUE, as.sp=TRUE)
}
\arguments{
  \item{X}{
    The point pattern.  An object of class \code{\link{sp2}} or
    \code{\link{sp3}}, or \code{\link{matrix}}. 
  }
  \item{win}{
    The window that contains the point pattern, must be an object of
    class \code{\link{sphwin}} with \code{type=="bandcomp"}.  Only used
    if \code{X} is an object of class \code{\link{matrix}}. 
  }
  \item{upper,lower}{
    Logicals, denoting whether the output should contain the points in
    the upper and lower caps respectively. 
  }
  \item{as.sp}{
    Optional.  If \code{TRUE}, and \code{X} is an object of class
    \code{\link{sp2}} or \code{\link{sp3}}, the output will be the an
    object or list of objects the same class as \code{X}. 
  }
}
\details{
  A complement of a band is the union of two disjoint caps with parallel
  bases.  This code separates the points in such a window by the cap
  that they are in.  In this context the upper cap is the cap with base
  defined by the first entry in the argument \code{param} in the
  \code{sphwin} object that defines the window, the lower cap with base
  defined by the second argument. 
}
\value{
  If \code{upper=TRUE} and \code{lower=FALSE} a matrix containing the
  locations of the points in the upper cap (if also \code{as.sp=FALSE})
  or an object of the same class as \code{X}. 

  If \code{upper=FALSE} and \code{lower=TRUE} a matrix containing the
  locations of the points in the lower cap (if also \code{as.sp=FALSE})
  or an object of the same class as \code{X}. 

  If \code{upper=TRUE} and \code{lower=TRUE} a list containing the
  locations of the points in the upper cap as matrices (if also
  \code{as.sp=FALSE}) or objects of the same class as \code{X}. 
}
\author{
Tom Lawrence <email:tjlawrence@bigpond.com>
}
\seealso{
  \code{\link{gcdist}}
}
\examples{
X <- rpoispp.sphwin(lambda=10, win=sphwin("bandcomp", param=c(pi/3, 2*pi/3)))
Xupper <- splitcaps(X=X, upper=TRUE, lower=FALSE)
Xlower <- splitcaps(X=X, upper=FALSE, lower=TRUE)
Xboth <- splitcaps(X=X, upper=TRUE, lower=TRUE)
}
\keyword{spatial}
\keyword{manip}

