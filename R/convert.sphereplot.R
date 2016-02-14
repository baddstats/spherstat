convert.sphereplot <- function(X, inverse=FALSE, sp.dim=2) {
	if(inverse) {
		stopifnot(sp.dim==2 || sp.dim==3)
		stopifnot(inherits(X, c("matrix", "numeric")))
		X <- cbind(90-X[,1], X[,2])
		X <- X*(180/pi)
		if(sp.dim==3) {
			X <- convert3(X)
		}
	}
	else {
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
		X <- X*(180/pi)
		X <- cbind(X[,2], (90-X[,1]))
	}
	X
}
