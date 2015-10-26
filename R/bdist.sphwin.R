bdist.sphwin <-
function(X, win=sphwin(type="sphere")) {
stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
if(inherits(X, "matrix")) {stopifnot(inherits(win, "sphwin"))} else {
win <- X$win
X <- X$X
stopifnot(inherits(win, "sphwin"))
}
n <- nrow(X)
if(nrow(X) > 0) {bdists <- switch(win$type,
sphere = {rep(Inf, n)},
band = {mindist.band(X=rot.sphere(points=rot.sphere(points=X, northpole=win$ref), northpole=c(0, win$param[2])) , win=win)},
bandcomp = {mindist.band(X=rot.sphere(points=rot.sphere(points=X, northpole=win$ref), northpole=c(0, win$param[2])) , win=win)},
wedge = {
sph.poly <- sphwin(type="polygon", param=matrix(c(0,0, pi/2, 0, pi, 0, pi/2, win$param[1], 0, 0), nrow=5, ncol=2, byrow=TRUE), ref = rep(0,4), rad=win$rad)
bd <- mindist.polygon(X=rot.sphere(points=X, northpole=win$ref, inverse=TRUE), win=sph.poly)
bd},
polygon = {bd <- mindist.polygon(X=X, win=win)
bd},
quadrangle = {Xdists <- cbind(X[,1]-win$param[1], win$param[2]-X[,1], X[,2], win$param[3]-X[,2])
bd <- apply(Xdists, 1, min)
bd
},
{stop("Unrecognised window type")})
} else {bdists <- NA}
bdists
}
