rot.sphere <- function(points, northpole, rad=1, inverse=FALSE) {
  if(inherits(points, c("sp2", "sp3"))) {
    points1 <- points$X
    rad <- points$win$rad
  } else if(is.matrix(points)) {
    points1 <- points
  } else if(is.numeric(points) && length(points) %in% c(2,3)) {
    points1 <- matrix(points, nrow=1)
  } else stop("points should be a point pattern or a matrix")
  if(ncol(points1)==2){
    points1 <- convert3(points1)
  }
  rotmat <- rot.matrix(northpole)
  if(inverse) {
    rotpoints <- t(rad*(t(rotmat) %*% t(points1)))
  } else {
    rotpoints <- t(rad*(rotmat %*% t(points1)))
  }
  output <- convert2(rotpoints)
  output
}
