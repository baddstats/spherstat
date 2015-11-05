#'  ptinsphpoly.R
#'
#' Test whether points lie in a spherical polygon
#'

ptinsphpoly <- function(X, win, P) {
  stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
  stopifnot(inherits(P, c("sp2", "sp3", "matrix")))
  stopifnot(inherits(win, "sphwin") && win$type == "polygon")
  rad <- win$rad
  if(inherits(X, "sp3") || (is.matrix(X) && ncol(X) == 3)) X <- convert2(X, rad)
  if(inherits(P, "sp3") || (is.matrix(P) && ncol(P) == 3)) P <- convert2(P, rad)
  if(inherits(X, "sp2")) X <- X$X
  if(inherits(P, "sp2")) P <- P$X
  nx <- nrow(X)
  Xlat <- pi/2 - X[,1]
  Xlon <- X[,2]
  Plat <- pi/2 - P[,1]
  Plon <- P[,2]
  vlat <- pi/2 - win$param[,1]
  vlon <- win$param[,2]
  # remove duplicated vertex
  nv <- nrow(win$param) - 1
  vlat <- vlat[1:nv]
  vlon <- vlon[1:nv]
  # go
  zz <- .C("RcallPtInSphPoly",
           plat = as.double(Xlat),
           plon = as.double(Xlon),
           np = as.integer(nx), 
           vlat = as.double(vlat),
           vlon = as.double(vlon),
           nv = as.integer(nv),
           xlat = as.double(Plat),
           xlon = as.double(Plon),
           plocation = as.integer(integer(nx)))
  (zz$plocation == 1)
}
