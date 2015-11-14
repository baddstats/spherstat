convert.globe <- function(X) {
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
	Xglobe <- X*(180/pi)
	Xglobe.out <- cbind(Xglobe[,2]-180, 90-Xglobe[,1])
	Xglobe.out
}
