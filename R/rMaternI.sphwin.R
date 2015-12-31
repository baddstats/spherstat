rMaternI.sphwin <- function(kappa, r, win=sphwin(type="sphere"), nsim=1, drop=TRUE, stationary=TRUE, as.sp=TRUE, ndim="2") {
	stopifnot(r > 0)
	stopifnot(kappa > 0)
	stopifnot(sum(ndim==c("2", "3")) == 1)
	output <- list()
	if(stationary==TRUE) {
		wina <- sphwin(type="sphere", rad=win$rad)
	}
	else {
		wina <- win
	}
	for(i in 1:nsim) {
		X1 <- rpoispp.sphwin(lambda=kappa, win=wina, as.sp=FALSE)
		nrX <- nrow(X1)
		X2 <- runif.sphwin(n=1, win=wina, as.sp=FALSE)
		for(j in 1:(nrX-1)) {
			Xtest <- runif.sphwin(n=1, win=wina, as.sp=FALSE)
			if(min(nncrosssph(Xtest, X2)) > r) {
				X2 <- rbind(X2, Xtest)
			}
		}
		if(stationary==TRUE) {
			Xin <- in.W(X2, win)
			X2 <- X2[Xin, , drop=FALSE]
		}
		output[[i]] <- switch(ndim,
			"2" = {
				X2
			},
			"3" = {
				convert3(X2)
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
