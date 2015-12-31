rMaternII.sphwin <- function(kappa, r, win=sphwin(type="sphere"), nsim=1, drop=TRUE, stationary=TRUE, as.sp=TRUE, ndim="2") {
	stopifnot(sum(ndim==c("2","3")) == 1)
	stopifnot(inherits(win, "sphwin"))
	stopifnot(kappa > 0)
	stopifnot(r > 0)
	if(stationary) {
		wina <- sphwin(type="sphere")
	}
	else {
		wina <- win
	}
	output <- list()
	for(i in 1:nsim) {
		X1 <- rpoispp.sphwin(lambda=kappa, win=wina, as.sp=FALSE)
		nrX <- nrow(X1)
		X2 <- runif.sphwin(n=1, win=wina, as.sp=FALSE)
		for(k in 1:(nrX-1)) {
			Xtest <- runif.sphwin(n=1, win=wina, as.sp=FALSE)
			to.keep <- c()
			nrX2 <- nrow(X2)
			if(length(X2)==0) {
				## cat("0 ")
				X2 <- rbind(X2, Xtest)
			}
			else {
				for(j in 1:nrX2) {
					if(gcdist(Xtest, matrix(X2[j,], nrow=1, ncol=2, byrow=TRUE)) > r) {
						to.keep <- c(to.keep, j)
					}
				}
				if(length(to.keep) == nrX2) {
					X2 <- rbind(X2, Xtest)
				}
				else {
					X2 <- X2[to.keep,,drop=FALSE]
				}
			}
		}
		if(stationary) {
			X2in <- in.W(X2, win)
			X2 <- X2[X2in, , drop=FALSE]
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
