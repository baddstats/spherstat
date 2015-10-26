pairdistsph <- function(X, rad=1) {
stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
if(inherits(X, "matrix")) {X <- X} else {
rad <- X$win$rad
X <- X$X}
nrX <- nrow(X)
if(ncol(X)==2) {X <- convert3(X)}
if(nrX==0) {d <- NA} else {
        Y <- X %*% t(X)
        d <- rad * acos(pmin(1, pmax(-1, Y/rad)))}
d <- matrix(d, nrow=nrX, ncol=nrX)    
return(d)
      }