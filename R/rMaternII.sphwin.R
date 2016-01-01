rMaternII.sphwin <- function(kappa, r, win=sphwin(type="sphere"), nsim=1, drop=TRUE, stationary=TRUE, as.sp=TRUE, ndim="2") {
	stopifnot(sum(ndim==c("2","3")) == 1)
	stopifnot(inherits(win, "sphwin"))
	stopifnot(kappa > 0)
	stopifnot(r > 0)
	wina <- if(stationary) sphwin(type="sphere") else win
	output <- list()
	for(i in 1:nsim) {
          X <- rpoispp.sphwin(lambda=kappa, win=wina, as.sp=FALSE)
          ## find all pairs of r-close points
          close <- (pairdistsph(X) < r)
          diag(close) <- FALSE
          ## attach arrival times to points.
          nrX <- nrow(X)
          tX <- runif(nrX)
          ## make a matrix which is TRUE if j arrived earlier than i
          earlier <- outer(tX, tX, ">")
          ## point i is killed if there is any earlier point j that is r-close
          kill <- apply(close & earlier, 1, any)
          ## delete killed points
          X <- X[!kill, , drop=FALSE]
          if(stationary) {
            Xin <- in.W(X, win)
            X <- X[Xin, , drop=FALSE]
          }
          output[[i]] <-
            switch(ndim,
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
	if(nsim==1 && drop) {
          output <- output[[1]]
	}
	output
}
