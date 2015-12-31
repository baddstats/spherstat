rpoispp.sphwin <- function(lambda, win=sphwin(type="sphere"), lmax=NULL, ..., nsim=1, drop=TRUE,  as.sp=TRUE, ndim="2") {
	stopifnot(inherits(win, "sphwin"))
	stopifnot(sum(ndim==c("2","3")) == 1)
	output <- list()
	for(i in 1:nsim) {
		if(is.numeric(lambda) && length(lambda)==1 && lambda>0) {
			asph <- area.sphwin(w=win)
			stopifnot(asph > 0)
			n <- rpois(1, lambda*asph)
			output[[i]] <- runif.sphwin(n, win=win, as.sp=as.sp, ndim=ndim)
		}
		else if (is.function(lambda) && !is.null(lmax)) {
	        	X <- rpoispp.sphwin(lambda=lmax, win=win, as.sp=as.sp, ndim=ndim)
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
			output[[i]] <- switch(ndim,
				"2" = {
					X
				},
				"3" = {
					convert3(X)
				},
				stop("ndim unrecognised.  See help file for details.")
			)
			if(as.sp) {
				if(ndim=="2") {
					output[[i]] <- sp2(output[[i]], win)
				}
				else {
					output[[i]] <- sp3(output[[i]], win)
				}
			}
		}
		else {
			stop("lambda and/or lmax improperly defined.  See help page for information on their how to define them.")
		}
	}
	if(nsim==1 && drop==TRUE) {
		output <- output[[1]]
	}
	output
}
