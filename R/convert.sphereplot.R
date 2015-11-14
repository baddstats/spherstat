convert.sphereplot <- function(X) {
	stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
	if(inherits(X, "matrix")) {
		X <- X
	}
	else {
	X <- X$X
	}
	if(ncol(X)==3) {
		X <- convert2(X)
	}
	Xsp <- X*(180/pi)
	Xsp.out <- cbind(Xsp[,2], (90-Xsp[,1]))
	Xsp.out
}
