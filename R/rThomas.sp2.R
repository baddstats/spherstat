rThomas.sp2 <- function(alpha, scale, mu, win=sphwin(type="sphere"), parents=FALSE) {
	stopifnot(inherits(win, "sphwin") && alpha > 0 && scale > 0 && mu > 0)
	rp <- rpoispp.sp2(lambda=alpha, win=win, as.sp=FALSE)
	rpl <- nrow(rp)
	rThom2 <- rp
	for(i in 1:rpl) {
		np <- rpois(1, mu)
		if(np > 0) {
			rThom1 <-  rFisher(n=np, mode=rp[i,], kappa=scale, win=sphwin(type="sphere", rad=win$rad))$X
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
	output <- sp2(X=output, win=win)
	output
}
