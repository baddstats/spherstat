rHardcore.sphwin <- function(beta, R=0, p=0.5, m=100, win=sphwin(type="sphere"), proper=TRUE, nsim=1, drop=TRUE, expand=TRUE, as.sp=TRUE, ndim="2") {
	stopifnot(beta > 0)
	stopifnot(p >= 0 && p <= 1)
	stopifnot (R >= 0)
	stopifnot(m > 0)
	stopifnot(inherits(win, "sphwin"))
	stopifnot(sum(ndim==c("2", "3"))==1)
	if(expand) {
		wina <- sphwin(type="sphere", rad=win$rad)
	}
	else {
		wina <- win
	}
	output <- list()
	for(i in 1:nsim) {
		X <- rpoispp.sphwin(lambda=beta, win=wina, as.sp=FALSE)  
		n <- nrow(X)
		if(proper) {
			while(sort(gcdist(X,X))[n+1] <= R) {
				prop <- runif(1,0,1)
				if(prop <= p) {
					X <- birth.sphwin(X=X, beta=beta, gamma=0, R=R, p=p, n=n, win=wina)
				}
				else {
					X <- death.sphwin(X=X, beta=beta, gamma=0, R=R, p=p, n=n, win=wina)
				}
				n <- nrow(X)
			}
		} else {
			for(i in 1:m) {
				prop <- runif(1,0,1)
				if(prop <= p) {
					X <- birth.sphwin(X=X, beta=beta, gamma=0, R=R, p=p, n=n, win=wina)
				}
				else {
					X <- death.sphwin(X=X, beta=beta, gamma=0, R=R, p=p, n=n, win=wina)
				}
				n <- nrow(X)
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
	if(drop==TRUE && nsim==1) {
		output <- output[[1]]
	}
	output
}

