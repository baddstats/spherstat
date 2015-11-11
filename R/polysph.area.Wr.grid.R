polysph.area.Wr.grid <- function(points, win, r=0, nlon=100) {
	stopifnot(inherits(points, c("sp2", "sp3", "matrix")))
	if(inherits(points, "matrix")) {
		stopifnot(inherits(win, "sphwin"))
	}
	else {
		win <- points$win
		points <- points$X
		stopifnot(inherits(win, "sphwin"))
	}
	posstypes <- c("sphere", "band", "wedge", "polygon", "quadrangle")
	if(!(win$type %in% posstypes)) {
		stop(paste("Type must be one of "), paste(posstypes, collapse=" or "))
	}
	rad <- win$rad
	lr <- length(r)
	lg <- nrow(points)
	rmat <- matrix(rep(r, lg), nrow=lr, ncol=lg, byrow=FALSE)
	cellarea <- (points[2,2]-points[1,2])*(rad^2)*(cos(points[1,1])-cos(points[(nlon+1),1]))
	if(win$type=="quadrangle") {
		distmat <- matrix(rep(bdist.sphwin(X=points, win=win), lr), nrow=lr, ncol=lg, byrow=TRUE)
		ndists <- rowSums(distmat >= rmat)
	}
	else {
		isinW <- in.W(points=points, win=win)
		numinW <- sum(isinW)
		isinW1 <- 1/as.numeric(isinW)
		isinWmat <- matrix(rep(isinW1), nrow=lr, ncol=lg, byrow=TRUE)
		distmat <- matrix(rep(bdist.sphwin(X=points, win=win), lr), nrow=lr, ncol=lg, byrow=TRUE)
		disttest <- distmat*isinWmat
		ndists <- rowSums(disttest >= rmat) - (lg-numinW)
	}
	areaWr <- ndists*cellarea
	areaWr
}