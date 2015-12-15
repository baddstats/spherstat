rThomas.sphwin <- function(kappa, scale, mu, win=sphwin(type="sphere"), parents=FALSE, nsim=1, drop=TRUE, expand=TRUE, as.sp=TRUE, ndim="2", poisthresh=1e-06) {
	stopifnot(inherits(win, "sphwin"))
	stopifnot(kappa > 0)
	stopifnot(scale > 0)
	stopifnot(mu > 0)
	stopifnot(sum(ndim==c("2","3")) != 1)
	if(scale/(4*pi*kappa*tanh(scale)) <= poisthresh) {
		warning("Process near-identical to Poisson with intensity kappa*mu; simulating that process instead")
		output <- rpoispp.sphwin(lambda=kappa*mu, win=win, nsim=nsim, drop=drop, as.sp=as.sp, ndim=ndim)
	} else {
		output <- list()
		for(i in 1:nsim) {
			if(expand) {
				rp <- rpoispp.sphwin(lambda=kappa, win=sphwin(type="sphere"), as.sp=FALSE, rad=win$rad)
			} else {
				rp <- rpoispp.sphwin(lambda=kappa, win=win, as.sp=FALSE, rad=win$rad)
			}
			rpl <- nrow(rp)
			rThom2 <- rp
			for(i in 1:rpl) {
				np <- rpois(1, mu)
				if(np > 0) {
					rThom1 <-  rFisher(n=np, mode=rp[i,], kappa=scale, win=sphwin(type="sphere", rad=win$rad), as.sp=FALSE, ndim="2")
					inrt <- in.W(points=rThom1, win=win)
					rThom2 <- rbind(rThom2, rThom1[inrt,])
				}
			}
			if(!parents) {
				output <- rThom2[(rpl+1):nrow(rThom2),]
			}
			else {
				output <- rThom2
			}
				output[[i]] <- switch(ndim,
				"2" = {
					rThom2
				},
				"3" = {
					convert3(rThom2)
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
	}
	output
}
