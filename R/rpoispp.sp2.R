rpoispp.sp2 <- function(lambda, win=sphwin(type="sphere"), lmax=NULL, as.sp=TRUE, sp.dim="2", ...) {
	stopifnot(inherits(win, "sphwin"))
	if(is.numeric(lambda) && length(lambda)==1 && lambda>0) {
		asph <- area.sphwin(w=win)
		stopifnot(asph > 0)
		n <- rpois(1, lambda*asph)
		X <- runif.sphwin(n, win=win, as.sp=as.sp, sp.dim=sp.dim)
	}
	else if (is.function(lambda) && !is.null(lmax)) {
	        X <- rpoispp.sp2(lambda=lmax, win=win, as.sp=as.sp, sp.dim=sp.dim)
		if(inherits(X, c("sp2", "sp3"))) {
			X <- X$X
		}
		nrX <- nrow(X)
        	if(nrX != 0) {
        		prob <- lambda(X, ...)/lmax
        		u <- runif(nrX)
        		retain <- (u <= prob)
        		X <- X[retain, ]
		}
		if(as.sp==TRUE) {
			X <- switch(sp.dim,
				"2" = {
					sp2(X, win)
				},
				"3" = {
					sp3(X, win)
				},
				stop("sp.dim incorrectly defined")
			)
		}
	}
	else {
		stop("lambda and/or lmax improperly defined.  See help page for information on their how to define them.")
	}
	X
}
