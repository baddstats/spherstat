crossdistsph <- 
function(X, Y, rad=1) {
stopifnot(inherits(X, c("sp2", "sp3", "matrix")) && inherits(Y, c("sp2", "sp3", "matrix")))
if(inherits(X, c("sp2", "sp3")) && inherits(Y, c("sp2", "sp3")) && X$win$rad != Y$win$rad) {stop("Objects have different rad values")}
if(inherits(X, "matrix")) {X <- X} else {
rad <- X$win$rad
X <- X$X
}
if(inherits(Y, "matrix")) {Y <- Y} else {
if(!inherits(Y, c("sp2", "sp3"))) {rad <- Y$rad}
Y <- Y$X}
nrX <- nrow(X)
nrY <- nrow(Y)
if(ncol(X)==2) {X <- convert3(X)}
if(ncol(Y)==2) {Y <- convert3(Y)}
if(nrX==0 && nrY==0) {d <- NA} else {
        XY <- X %*% t(Y)
        d <- rad * acos(pmin(1, pmax(-1, XY/rad))/rad)}
d <- matrix(d, nrow=nrX, ncol=nrY)         
return(d)
      }
