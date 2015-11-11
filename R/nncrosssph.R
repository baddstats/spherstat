nncrosssph <- function(X, Y, rad=1) {
	stopifnot(inherits(X, c("sp2", "sp3", "matrix")) && inherits(Y, c("sp2", "sp3", "matrix")))
	if(inherits(X, "matrix")) {
		X<-X
	}
	else {
		rad <- X$win$rad
		X <- X$X
	}
	if(inherits(Y, "matrix")) {
		Y<-Y
	}
	else {
		if(Y$win$rad != rad) {
			stop("X and Y have different radii")
		}
	Y <- Y$X
	}
	sphdist <- crossdistsph(X=X, Y=Y)
	sphmat <- apply(sphdist, 2, min)
	sphmat
}