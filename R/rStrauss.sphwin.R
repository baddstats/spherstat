rStrauss.sphwin <- function(beta, gamma, R=0, p=0.5, m=100, win=sphwin(type="sphere"), nsim=1, drop=TRUE, expand=TRUE, as.sp=TRUE, ndim=2) {
	stopifnot(inherits(win, "sphwin"))
	stopifnot(beta > 0)
	stopifnot(gamma >=0 && gamma <= 1)
	stopifnot(R >= 0)
	stopifnot(p >=0 && p <= 1)
	stopifnot(m > 0)
	stopifnot(sum(ndim==c("2","3")) != 1)
	if(gamma==0) {
		cat("Warning: Since gamma=0, simulated pattern from a Hard-Core process, not a Strauss process.")
	}
	if(expand) {
		wina <- sphwin(type="sphere", rad=win$rad)
	}
	else {
		wina <- win
	}
	output <- list()
	for(i in 1:nsim) {
		X <- rpoispp.sphwin(beta, win=wina, as.sp=FALSE)
		for (i in 1:m) {
			n <- nrow(X)
		  	prop <- runif(1,0,1)
		 	if(prop <= p) {
		  		X <- birth.sphwin(X=X, beta=beta, gamma=gamma, R=R, p=p, n=n, win=wina)
			} else {
		  		X <- death.sphwin(X=X, beta=beta, gamma=gamma, R=R, p=p, n=n, win=wina)
			}
		}
		if(expand) {
			inX <- in.W(X, win)
			X <- X[inX,]
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
	if(nsim==1 && drop==TRUE) {
		output <- output[[1]]
	}
	output
}
