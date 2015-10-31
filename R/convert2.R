convert2 <- 
function(points, rad=1) {
  # ajb: in the next line, 'X' has not yet been defined
stopifnot(length(points) ==3 || ncol(points) == 3 || inherits(X, "sp3"))
X <- points
if(inherits(X, "sp3")) {
rad <- points$win$rad
points <- points$X
}
if(!inherits(points, "matrix")) {points1 <- matrix(points, ncol=3)/rad} else {points1 <- points/rad}
# ajb: 'n' is defined but not used 
n <- nrow(points1)
# ajb: is it correct to divide by 'rad' here again?
# ajb: The expression points1[,3]/rad is re-evaluated several times.
theta <- ifelse(1-(points1[,3]/rad) <= 0,  0, ifelse(1+(points1[,3]/rad) <= 0, pi, acos(cround(points1[,3]/rad))))
stheta <- sin(theta)
phi <- ifelse(abs(stheta)<10^-15, 0, atan2((points1[,2]/stheta), (points1[,1]/stheta))%%(2*pi))
output <- cbind(theta, phi)
if(inherits(X, "sp3")) {output <- sp2(X=output, win=X$win)}
output
}
