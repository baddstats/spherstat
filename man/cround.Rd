\name{triground}
\alias{cround}
\alias{sround}
\title{
  Rounding for Trigonometric Functions
}
\description{
  Round the input for or output from trigonometric functions, to
  overcome some issues relating to rounding error. 
}
\usage{
cround(x)
sround(x)
}
\arguments{
  \item{x}{
    The number or numbers to be rounded.  Can be a vector, matrix, numeric, list or data frame.
  }
}
\details{
  The function \code{cround} rounds a number such that it is between -1
  and 1.  This overcomes rounding error whereby a number that has
  absolute value of 1 or very slightly less, is calculated by R to have
  absolute value of very slightly larger than 1, which can trigger
  errors if that calculated value is automatically passed to
  e.g. \code{\link{acos}} or \code{\link{asin}}.  The function
  \code{sround} rounds a number that is extremely close to 0, to 0; this
  overcomes a minor calculation error in some R functions
  (e.g. trigonometric functions) and thus may trigger errors in some R
  functions (e.g. \code{\link{acos}} or \code{\link{asin}}, if the value
  being rounded is to be added to or subtracted from 1 or -1. 
}
\value{
  The rounded values, in the same format as the input.
}

\author{
  Tom Lawrence <email: tjlawrence@bigpond.com>
}

\section{Warning }{
  \code{cround} will round any number greater than 1 to
  1; and any number less than -1 to -1.  Care should be taken when using
  \code{cround} to ensure that its input are between (or close to one of) -1
  and 1, that is that genuinely erroneous values are identified as such.  
}
\examples{
  x <- c(-1-(10^17), -1+(10^-17), -0.5, 0, 0.5, 1-(10^-17), 1+(10^-17))
  cround(x)

  y <- c(-0.5, -10^-17, 10^-17, 0.5)
  sround(y)
}
\keyword{arith}

