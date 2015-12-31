rFisher <- function(n, mode, kappa, win=sphwin(type="sphere"), nsim=1, drop=TRUE, as.sp=TRUE, ndim="2") {
	stopifnot(inherits(win, "sphwin"))
	stopifnot(kappa > 0)
	stopifnot(sum(ndim==c("2", "3"))==1)
	stopifnot(n > 0)
	lambda <- exp(-2*kappa)
	output <- list()
	for(i in 1:nsim) {
		r1 <- runif(n,0,1)
		r2 <- runif(n,0,1)
		theta <- 2*asin(sqrt(-log(lambda + (1-lambda)*r1)/(2*kappa)))
		theta <- theta
		phi <- 2*pi*r2
		rf1 <- cbind(theta, phi)
		rf1 <- rot.sphere(points=rf1, northpole=mode, inverse=TRUE)
		rf1in <- in.W(rf1, win)
		r2keep <- rf1[rf1in,, drop=FALSE]
		output[[i]] <- switch(ndim,
			"2" = {
				r2keep
			},
			"3" = {
				convert3(r2keep)
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
	if(nsim==1 && drop==TRUE) {
		output <- output[[1]]
	}
	output
}
