pairdistsph <- function(X, rad=1) {
  stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
  if(!inherits(X, "matrix")) {
    rad <- X$win$rad
    X <- X$X
  }
  nrX <- nrow(X)
  if(nrX == 0) return(matrix(, nrow=0, ncol=0))
  if(ncol(X)==2) {X <- convert3(X, rad)}
  Y <- X %*% t(X)
  if(rad == 1) {
    d <- acos(pmin(1, pmax(-1, Y))) }
  else {
    d <- rad * acos(pmin(1, pmax(-1, Y/rad)))
  }
  d <- matrix(d, nrow=nrX, ncol=nrX)    
  return(d)
}
