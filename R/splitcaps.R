splitcaps <- function(X, win=sphwin(type="sphere"), upper=TRUE, lower=TRUE, as.sp=TRUE) {
	stopifnot(inherits(X, c("sp2", "sp3", "matrix")) && (upper+lower)!=0)
	X0 <- X
	if(inherits(X, "matrix")) {stopifnot(inherits(win, "sphwin") && win$type=="bandcomp")} else {
		win <- X$win
		X <- X$X
		stopifnot(inherits(win, "sphwin") && win$type=="bandcomp")
	}
	Xrefdist1 <- acos(cround(convert3(X) %*% t(convert3(win$ref))))
	Xrefdist2 <- acos(cround(convert3(X) %*% t(convert3(c(pi-win$ref[1], (win$ref[2]+pi)%%(2*pi))))))
	if(upper) {
		Xupper <- X[Xrefdist1 <= win$param[1]*win$rad,]
		if(as.sp) {
			if(inherits(X0, "sp2")) {
				Xupper <- sp2(Xupper, win=sphwin(type="band", param=c(0, win$param[1]), ref=win$ref, rad=win$rad))
			} else {
				Xupper <- sp3(Xupper, win=sphwin(type="band", param=c(0, win$param[1]), ref=win$ref, rad=win$rad))
			}
		}
	}
	if(lower) {
		Xlower <- X[Xrefdist2 <= (pi-2*win$param[1])*win$rad,]
		if(as.sp) {
			if(inherits(X0, "sp2")) {
				Xlower <- sp2(Xlower, win=sphwin(type="band", param=c(win$param[2], pi), ref=win$ref, rad=win$rad))
			} else {
				Xlower <- sp3(Xlower, win=sphwin(type="band", param=c(win$param[2], pi), ref=win$ref, rad=win$rad))
			}
		}
	}
	if(upper) {
		if(lower) {
			Xout <- list(Xupper=Xupper, Xlower=Xlower)
		} else {
			Xout <- Xupper
		}
	} else if(lower) {
		Xout <- Xlower
	}
	Xout
}