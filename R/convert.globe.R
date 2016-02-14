convert.globe <- function(X, inverse=FALSE, sp.dim=2) {
	if(inverse) {
		stopifnot(sp.dim==2 || sp.dim==3)
		stopifnot(inherits(X, c("list", "matrix"))  & !inherits(X, c("sp2", "sp3")))
		if(inherits(X, "list")) {
			Xlon <- X[[1]]
			Xlat <- X[[2]]
		}
		else {
			Xlon <- X[,1]
			Xlat <- X[,2]
		}
		Xlon <- ifelse(Xlon >= 0, Xlon, Xlon+360)
		X <- cbind(90-Xlat, Xlon)*(pi/180)
		if(sp.dim==3) {
			X <- convert3(X)
		}
	}
	else {
		stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
		if(!inherits(X, "matrix")) {
			X <- X$X
		}
		if(ncol(X)==3) {
			X <- convert2(X)
		}
		Xdegrees <- X*(180/pi)
		Xlat <- 90 - Xdegrees[,1]
		Xlon <- Xdegrees[,2]
		Xlon <- ifelse(Xlon <= 180, Xlon, Xlon - 360)
		X <- list(lon=Xlon, lat=Xlat)
	}
	return(X)
}
