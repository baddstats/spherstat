\name{galaxies}
\alias{galaxies}
\docType{data}
\title{
  Galaxy Survey
}
\description{
  The sky positions of 10610 galaxies from the New General Catalogue.
}
\usage{data("galaxies")}
\format{
  A numeric matrix with two columns containing the longitude (in degrees
  from 0 to 360) and latitude (in degrees from -90 to 90).
}
\details{
  From the original list of 13957 galaxies given by Steinicke (2015)
  we have retained only those
  which lie in the region of space that is not partially
  obscured by the Milky Way. The obscured region is taken to be the
  region lying between 
}
\source{
  Steinicke (2015)
}
\references{
  Steinicke, W. (2015) Revised New General Catalogue and Index
  Catalogue.
  Available from \url{www.klima-luft.de/steinicke/index_e.html}.
  Accessed 14 May 2015.
}
\examples{
data(galaxies)
}
\keyword{datasets}
