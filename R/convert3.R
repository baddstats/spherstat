convert3 <-
function(points, rad=1) {
stopifnot(length(points) == 2 || ncol(points)==2 || inherits(points, "sp2"))
points1 <- points
# ajb: 'points1' is defined but never used
if(inherits(points, "sp2")) {
rad <- points$win$rad
points <- points$X
}
# ajb: next line coul generate an error: class(points) could be a vector
if(class(points) != "matrix") {points <- matrix(points, ncol=2, byrow=TRUE)} else {points <- points}
stopifnot(ncol(points)==2)
# ajb: 'n' is defined but not used
n <- nrow(points)
# ajb: sin(points[,1]) is re-evaluated several times
points3 <-  cbind(sin(points[,1]) * cos(points[,2]),
                       sin(points[,1]) * sin(points[,2]),
                       cos(points[,1]))
points3 <- sround(rad * points3)
points3
}
