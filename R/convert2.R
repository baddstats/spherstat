convert2 <- function(points, rad=1) {
	stopifnot(length(points) ==3 || ncol(points) == 3 || inherits(points, "sp3"))
	X <- points
	if(inherits(X, "sp3")) {
		rad <- points$win$rad
		points <- points$X
  	}
	if(!is.matrix(points)) {
		points <- matrix(points, ncol=3)/rad
	}
	theta <- acos(cround(points[,3]/rad))
	phi <- atan2(points[,2], points[,1]) %% (2*pi)
	output <- cbind(theta, phi)
	if(inherits(X, "sp3")) {
		output <- sp2(X=output, win=X$win)
	}
	output
}
