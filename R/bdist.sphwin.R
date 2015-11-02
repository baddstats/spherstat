bdist.sphwin <- function(X, win=sphwin(type="sphere")) {
  stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
  if(inherits(X, "matrix")) {stopifnot(inherits(win, "sphwin"))} else {
    win <- X$win
    X <- X$X
    if(!inherits(win, "sphwin"))
      stop("X$win should be a window")
  }
  n <- nrow(X)
  if(n == 0) return(numeric(0))
  rad <- win$rad
  radpar <- rad * win$param
  switch(win$type,
         sphere = {
           bdists <- rep(Inf, n)
         },
         bandcomp =,
         band = {
           centre <- matrix(win$ref, nrow=1)
           dcentre <- as.vector(gcdist(X, centre, rad=rad))
           bdists <- pmin(abs(dcentre - radpar[1]), abs(dcentre - radpar[2]))
         },
         wedge = {
           mat <- matrix(c(0,0, pi/2, 0, pi, 0, pi/2, win$param[1], 0, 0),
                         nrow=5, ncol=2, byrow=TRUE)
           sph.poly <- sphwin(type="polygon", param=mat,
                              ref = rep(0,4), rad=rad)
           Xrot <- rot.sphere(points=X, northpole=win$ref, inverse=TRUE)
           bdists <- mindist.polygon(Xrot, win=sph.poly)
         },
         polygon = {
           bdists <- mindist.polygon(X=X, win=win)
         },
         quadrangle = {
           Xdists <- cbind(X[,1]-win$param[1], win$param[2]-X[,1],
                           X[,2], win$param[3]-X[,2])
           bdists <- apply(Xdists, 1, min)
         },
         {stop("Unrecognised window type")})
  return(bdists)
}
