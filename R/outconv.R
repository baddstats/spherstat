outconv <- function(X, win, as.sp, sp.dim) {
	if(sp.dim=="3") {
		X <- convert3(X)
		if(as.sp==TRUE) {
			X <- sp3(X, win)
		}
	}
	else if(as.sp==TRUE) {
		X <- sp2(X=X, win=win)
	}
	X
}