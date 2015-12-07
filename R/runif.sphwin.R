runif.sphwin <- function(n, win=sphwin(type="sphere"), nsim=1, drop=TRUE, as.sp=TRUE, ndim="2") {
	stopifnot(inherits(win, "sphwin"))
	stopifnot(sum(ndim==c("2","3")) == 1)
	output <- list()
	for(i in 1:nsim) {
		X <- switch(win$type,
			sphere={
				runif.sphere(n=n, win=win)
			},
			band={
				runif.band(n=n, win=win)
			},
			bandcomp={
				runif.bandcomp(n=n, win=win)
			},
			wedge={
				runif.wedge(n=n, win=win)
			},
			polygon={
				runif.polygon(n=n, win=win)
			},
			quadrangle={
				runif.quadrangle(n=n, win=win)
			},
			stop("Unsupported shape type!")
		)
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
	if(nsim==1 && drop==TRUE) {
		output <- output[[i]]
	}
	output
}