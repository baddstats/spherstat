rMatClust.sphwin <- function(kappa, scale, mu, win=sphwin(type="sphere"), parents=FALSE, nsim=1, drop=TRUE,expand=TRUE, as.sp=TRUE, ndim="2") {
	stopifnot(inherits(win, "sphwin"))
	stopifnot(kappa > 0)
	stopifnot(mu > 0)
	stopifnot(scale > 0)
	stopifnot(sum(ndim==c("2", "3"))==1)
	rad <- win$rad
	output <- list()
	for(i in 1:nsim) {
		if(expand) {
			rp <- rpoispp.sphwin(lambda=kappa, win=sphwin(type="sphere", rad=win$rad), as.sp=FALSE)
		} else {
			rp <- rpoispp.sphwin(lambda=kappa, win=win, as.sp=FALSE)
		}
		rpl <- nrow(rp)
		rM <- rp
		for(j in 1:rpl) {
			nlam <- rpois(1, mu)
			if(nlam > 0) {
				daughtwin <- sphwin(type="band", param=c(0, scale/rad), ref=rp[i,])
				rMat1 <-  runif.sphwin(n=nlam,  win=daughtwin, as.sp=FALSE)
				inrm <- in.W(points=rMat1, win=win)
				rM <- rbind(rM, rMat1[inrm,,drop=FALSE])
			}
		}
		if(!parents) {
			X <- rM[(rpl+1):nrow(rM),,drop=FALSE]
		}
		else {
			X <- rM
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
	if(nsim==1 & drop==TRUE) {
		output <- output[[1]]
	}
	output
}
