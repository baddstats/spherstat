nndistsph <- function(X, rad=1) {
  stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
  if(!is.matrix(X)) {
    rad <- X$win$rad
    X <- X$X
  }
  if(nrow(X) == 0) return(numeric(0))
  if(ncol(X)==3) { X <- convert2(X, rad) }
  sphdist <- pairdistsph(X, rad)
  diag(sphdist) <- Inf
  sphmat <- apply(sphdist, 1, min)
  return(sphmat)
}
