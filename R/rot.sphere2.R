rot.sphere2 <- function(points, v, theta) {
	if(inherits(points, c("sp2", "sp3"))) {
		points1 <- points$X
	} else if(is.matrix(points)) {
		points1 <- points
	} else if(is.numeric(points) && length(points) %in% c(2,3)) {
		points1 <- matrix(points, nrow=1)
	} else stop("points should be a point pattern or a matrix")
	if(ncol(points1)==2){
		points1 <- convert3(points1)
	}
	if(inherits(v, c("sp2", "sp3"))) {
		v1 <- v$X
	} else if(is.matrix(v)) {
		v1 <- v
	} else if(is.numeric(v) && length(v) %in% c(2,3)) {
		v1 <- matrix(v, nrow=1)
	} else stop("v should be a point pattern with one point, a matrix with one row, or a numeric of length 2 or 3")
	if(ncol(v1)==2){
		v1 <- convert3(v1)
	}
	v1 <- v1/sqrt(sum(v1^2))
	ctheta <- cos(theta)
	stheta <- sin(theta)
	if(inherits(points1, "matrix")) {
		nrp <- nrow(points1)
		y <- matrix(nrow=nrp, ncol=3)
		for(i in 1:nrp) {
			y[i,] <- points1[i,]*ctheta + cross(v1,points1[i,])*stheta+v1*dot(v1,points1[i,])*(1-ctheta)
		}
	} else {		
		y <- points1*ctheta + cross(v1, points1)*stheta + v1*dot(v1,points1)*(1-ctheta)
	}
	if(inherits(points, "sp2")) {
		yout <- points
		yout$X <- convert2(y)
	} else if(inherits(points, "sp3")) {
		yout <- points
		yout$X <- y
	} else if(inherits(points, "matrix")) {
		ncp1 <- ncol(points)
		if(inherits(y, "numeric")) {
			y <- matrix(y, ncol=3, nrow=1)
		}
		if(ncp1==2) {
			yout <- convert2(cround(sround(y)))
		} else {
			yout <- cround(sround(y))
		}
	}
	else if(inherits(points, "numeric")) {
		yout <- y
		if(length(points)==2) {
			yout <- convert2(y)
		}
	}
	yout
}
