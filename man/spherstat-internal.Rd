\name{spherstat-internal}
\alias{birth.sphwin}
\alias{death.sphwin}
\alias{gridmat.nlon}
\alias{in.W.poly}
\alias{mindist.band}
\alias{mindist.polygon}
\alias{polysph.area.W}
\alias{polysph.area.Wr}
\alias{polysph.area.Wr.exact}
\alias{polysph.area.Wr.grid}
\alias{polysph.area.Wr.int}
\alias{runif.band}
\alias{runif.bandcomp}
\alias{runif.polygon}
\alias{runif.quadrangle}
\alias{runif.sphere}
\alias{runif.wedge}
% \alias{sphsincos.vert} DOES NOT EXIST
\alias{intcircs}
\alias{Kiso}
\alias{Kisoband}
\alias{Kisobc}
\alias{Kisocap}
\alias{Kisoengine}
\alias{Kisopoly}
\alias{Kisopolyengine}
\alias{Kisoquad}
\alias{Kisowedge}
\alias{Kisowedgeengine}
\alias{quadtopoly}
\alias{neighbourcount}
\alias{cumtab}
\alias{rmax.rule.sphwin}
\alias{ptinsphpoly}
\alias{allcoords}
\title{
Internal Spherstat Functions
}
\usage{
birth.sphwin(X, beta, gamma, R, p, n, rad, win)
death.sphwin(X, beta, gamma, R, p, n, win)
gridmat.nlon(colats, lons, ncolat, nlon)
in.W.poly(points, win)
mindist.band(X, win)
mindist.polygon(X, win)
polysph.area.W(win)
polysph.area.Wr(method, \dots)
polysph.area.Wr.exact(win, r)
polysph.area.Wr.grid(points, win, r, nlon)
polysph.area.Wr.int(win, r)
runif.band(n, win)
runif.bandcomp(n, win)
runif.polygon(n, win)
runif.quadrangle(n, win)
runif.sphere(n, win)
runif.wedge(n, win)
intcircs(n1, n2, h1, h2, tol)
Kiso(X, win, r, rad, Dmat, nrX, denom, lambda)
Kisoband(X, win, r, nrX, Dmat, disc, rad, denom, useC, lambda)
Kisobc(X, win, r, nrX, Dmat, disc, verbose, rad, denom, useC, lambda)
Kisocap(X, win, r, nrX, Dmat, disc, rad, denom, useC, lambda)
Kisoengine(xi3, xj3, win, ints, verbose, cdij)
Kisopoly(X, win, r, nrX, Dmat, disc, quadwin, rad)
Kisopolyengine(xi3, xj3, vk03, vk13, win, ints, cdij, nr, cvn, cv2n)
Kisoquad(X, win, r, nrX, Dmat, disc)
Kisowedge(X, win, r, nrX, Dmat, disc, verbose, rad)
Kisowedgeengine(xi, xj, ints, nr, c1, c2, ref, win, verbose)
quadtopoly(X, win)
neighbourcount(Dmat, r)
cumtab(z, m)
rmax.rule.sphwin(win)
ptinsphpoly(X, win, P)
allcoords(X)
}
\description{
Internal spherstat functions
}
\details{
These are generally not to be called by the user.
}
\author{
Adrian Baddeley and Tom Lawrence
}
\keyword{internal}
