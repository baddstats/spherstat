convert3 <-
function(points, rad=1) {
stopifnot(length(points) == 2 || ncol(points)==2 || inherits(points, "sp2"))
points1 <- points
if(inherits(points, "sp2")) {
rad <- points$win$rad
points <- points$X
}
if(class(points) != "matrix") {points <- matrix(points, ncol=2, byrow=TRUE)} else {points <- points}
stopifnot(ncol(points)==2)
n <- nrow(points)
points3 <-  cbind(sin(points[,1]) * cos(points[,2]),
                       sin(points[,1]) * sin(points[,2]),
                       cos(points[,1]))
points3 <- sround(rad * points3)
points3
}
