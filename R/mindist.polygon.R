mindist.polygon <-
function(X, win) {
stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
if(inherits(X, "matrix")) {stopifnot(inherits(win, "sphwin"))} else {
win <- X$win
X <- X$X
}
rad <- win$rad
n <- nrow(X)
lp <- nrow(win$param)
p3 <- convert3(win$param)
if(ncol(X) !=3) {X <- convert3(X)}
gc.pv1 <- gcdist(x=X, y=p3[1:lp,], rad=rad)
gc.pv2 <- t(gc.pv1[,1:(lp-1)])
gc.pv3 <- t(gc.pv1[,2:lp])
gc.v <- diag(gcdist(x=p3[1:(lp-1),], y=p3[2:lp,], rad=rad))
angs <- sphcos(d1=gc.pv2, d2=gc.v, d3=gc.pv3, theta=NULL, rad=rad)
dists <- sphsin(d1=gc.pv2, d2=NULL, theta1=pi/2, theta2=angs, rad=rad)
if(n==1) {md <- min(dists)} else {md <- apply(dists, 2, min)}
md
}
