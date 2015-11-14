death.sp2 <- function(X, beta, gamma, R, p, n, win) {
	rad <- win$rad
	xnum <- rdiscunif(1,1,n)
	if(xnum==1) {
		xvec <- 2:n
	}
	else if(xnum==n) {
		xvec <- 1:(n-1)
	}
	else {
		xvec=c(1:(xnum-1),(xnum+1):n)
	}
	X1 <- X[xvec,]
	xpoint <- matrix(X[xnum,], ncol=2, nrow=1, byrow=TRUE)
	xs <- ifelse(gcdist(xpoint, X1) <= R, 1, 0)
	deathcheck <- beta*(gamma^(sum(xs)-1))*(n/(1-p))*(p/(4*pi*rad^2))
	probdeath <- min(1, deathcheck)
	y <- runif(1,0,1)
	if(y >= probdeath) {
		Xfinal <- X
	}
	else {
		Xfinal <- X1
	}
	Xfinal
}
