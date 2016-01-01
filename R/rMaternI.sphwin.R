rMaternI.sphwin <- function(kappa, r, win=sphwin(type="sphere"), nsim=1, drop=TRUE, stationary=TRUE, as.sp=TRUE, ndim="2") {
	stopifnot(r > 0)
	stopifnot(kappa > 0)
	stopifnot(sum(ndim==c("2", "3")) == 1)
	output <- list()
	if(stationary) {
          wina <- sphwin(type="sphere", rad=win$rad)
	} else {
          wina <- win
	}
	for(i in 1:nsim) {
          X1 <- rpoispp.sphwin(lambda=kappa, win=wina, as.sp=FALSE)
          nnd <- nndistsph(X1)
          retain <- (nnd >= r)
          X2 <- X1[retain, , drop=FALSE]
          if(stationary) {
            Xin <- in.W(X2, win)
            X2 <- X2[Xin, , drop=FALSE]
          }
          output[[i]] <-
            switch(ndim,
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
	if(nsim==1 && drop) {
          output <- output[[1]]
	}
	output 
}
