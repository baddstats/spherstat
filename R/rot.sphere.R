rot.sphere <-
function(points, northpole, rad=1, inverse=FALSE) {
stopifnot(inherits(points, c("sp2", "sp3", "matrix")))
if(inherits(points, c("sp2", "sp3"))) {points1 <- points$X
rad <- points$win$rad} else {points1 <- points}
if(ncol(points1)==2){points1 <- convert3(points1)}
rotmat <- rot.matrix(northpole)
if(inverse) {rotpoints <- t(rad*(t(rotmat) %*% t(points1)))} else {rotpoints <- t(rad*(rotmat %*% t(points1)))}
output <- convert2(sround(rotpoints))
output
}
