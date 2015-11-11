sphwin <- function(type="sphere", param, ref=c(0,0), ref2=NULL, rad=1) {
	posstypes <- c("sphere", "band", "bandcomp", "wedge", "polygon", "quadrangle")
	if(!(type %in% posstypes)) {
		stop(paste("Type must be one of "), paste(posstypes, collapse=" or "))
	}
	stopifnot(is.finite(rad) && rad > 0)
	switch(type,
		sphere = {
			param=c()
			},
		band = {
			c(stopifnot(length(param)==2 && param[1] <= param[2] && length(ref)==2))
                },
		wedge = {
			c(stopifnot(length(param)==2 && length(ref)==2))
                },
	    	bandcomp = {
			c(stopifnot(length(param)==2 && param[1] <= param[2] && length(ref)==2))
                },
	    	polygon = {
			nrp <- nrow(param)
			c(stopifnot(ncol(param)==2 && param[1,]==param[nrp,]))
			if(!is.null(ref)) {
				c(stopifnot(length(ref)==nrp-1 && ref*(1-ref)==rep(0, nrp-1)))
			}
			else {
				ref <- rep(0, nrp-1)
			}
			if(!is.null(ref2)) {
				c(stopifnot(length(ref2)==length(ref) && !(any(ref2 < 0) || any(ref2 > pi/2))))
			}
			else {
				ref2 <- rep(pi/2, times=nrp-1)
			}
			p3 <- convert3(param)
			if(sum(cround(diag((p3[1:nrp,] %*% t(p3[1:nrp,])/(rad^2))[2:nrp,]))==-rad)>0){
				stop("Two consecutive vertices are diametrically opposite, hence the arc between them cannot be determined.  Please add a vertex between those vertices.")
			}
			if(sum(cround(diag((p3[1:nrp,] %*% t(p3[1:nrp,])/(rad^2))[2:nrp,]))==rad)>0){
				stop("Two consecutive vertices are at the same location.  One of these vertices must be removed in order for sphwin to create a polygonal window")
			}
                    },
		quadrangle = {
			stopifnot(length(param)==4 && param[1] <= param[2] && length(ref)==2)
		},
		{
			stop("Unsupported shape type!")
		}
       )
	if(type=="band" || type=="wedge") {
		if(ref[1]*(pi-ref[1])==0 && ref[2] != 0) {
			ref[2] <- 0
		}
	}
	result <- list(type=type, param=param, ref=ref, ref2=ref2, rad=rad)
	class(result) <- c("sphwin", class(result))
	result
}