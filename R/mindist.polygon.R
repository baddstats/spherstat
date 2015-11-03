mindist.polygon <- function(X, win) {
  stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
  if(inherits(X, "matrix")) {
    stopifnot(inherits(win, "sphwin"))
  } else {
    win <- X$win
    X <- X$X
  }
  rad <- win$rad
  nX <- nrow(X)
  if(nX == 0) return(Inf)
  if(ncol(X) !=3) {X <- convert3(X)}
  verts <- convert3(win$param)
  nvert <- nrow(verts)  ## note: last vertex = first vertex
  ## distance from each point to each vertex
  dXV <- gcdist(x=X, y=verts, rad=rad)
  dXVj      <- t(dXV[,-nvert,drop=FALSE])
  dXVjplus1 <- t(dXV[,-1,    drop=FALSE])
  ## polygon edge lengths
  elen <- gcdistPaired(x=verts[-nvert,,drop=FALSE],
                       y=verts[-1,,drop=FALSE],
                       rad=rad)
  ## angle X[i], V[j], V[j+1]
  ang1 <- sphcos(d1=dXVj, d2=elen, d3=dXVjplus1, theta=NULL, rad=rad)
  ## angle X[i], V[j+1], V[j]
  ang2 <- sphcos(d1=dXVjplus1, d2=elen, d3=dXVj, theta=NULL, rad=rad)
  ## perpendicular distance from X[i] to great circle containing edge
  perps <- sphsin(d1=dXVj, d2=NULL, theta1=pi/2, theta2=ang1, rad=rad)
  ## distance from X[i] to edge V[j], V[j+1]
  ##     If both angles are acute, use the perpendicular distance,
  ##     otherwise the smaller of the point-to-vertex distances
  dists <- ifelse(ang1 < pi/2 & ang2 < pi/2, perps, pmin(dXVj, dXVjplus1))
  md <- if(!is.matrix(dists)) min(dists) else apply(dists, 2, min)
  return(md)
}

