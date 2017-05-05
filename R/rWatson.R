## Function: rWatson
## Role: simulate datasets of n points from a Watson distribution with location parameter mode and concentration parameter kappa.   This distribution is Watson bipolar (kappa > 0), Watson girdle (kappa < 0), or uniform (kappa=0); this function calls the appropriate function for simulating that model
## Output: An object containing the simulated points; a list of datasets (unless drop=TRUE and nsim=1).  Each dataset is a matrix with 2 or 3 columns, or an object of class sp2 or sp3
## Arguments:
## n : number of points to be simulated (numeric of length 1, > 0)
## mode : location parameter, gives the location of a mode (Watson bipolar) or antimode (Watson girdle).  Used to determine which simulation function to call.
## kappa : concentration parameter, numeric of length 1
## win : the window in which to simulate the points
## squeeze : logical.  If TRUE, then uses a squeeze function to accelerate the simulation process.  Recommended only when kappa is large.
## nsim : number of datasets to simulate.  numeric of length 1
## drop : logical.  If TRUE and nsim=1, the output will not be a list.  Ignored if nsim !=1
## as.sp : logical.  If TRUE, the output is a (list of) sp2 or sp3 objects, depending on the value of ndim
## ndim : numeric, equal to 2 or 3.  This determines the locations of points should be given in 2 dimensions (i.e. polar coordinates) or 3 dimensions (i.e. Cartesian coordinates)

rWatson <- function(n, mode, kappa, win=sphwin(type="sphere"), squeeze=FALSE, nsim=1, drop=TRUE, as.sp=TRUE, ndim="2") {
	stopifnot(is.numeric(kappa))
	if(kappa > 0) {
		x <- rWatson.bipolar(n=n, mode=mode, kappa=kappa, win=win, squeeze=squeeze, nsim=nsim, drop=drop, as.sp=as.sp, ndim=ndim)
	} else if(kappa < 0) {
		x <- rWatson.girdle(n=n, mode=mode, kappa=kappa, win=win, squeeze=squeeze, nsim=nsim, drop=drop, as.sp=as.sp, ndim=ndim)
	} else {
		print("Note: Since kappa=0, this distribution is a uniform distribution and hence kappa, mode and squeeze a are ignored")
		x <- runif.sphwin(n, win=win, nsim=nsim, drop=drop, as.sp=as.sp, ndim=ndim)
	}
	x
}

## Function: rWatson.bipolar
## Role: simulate datasets of n points from a Watson bipolar distribution with location parameter mode.
## Output: An object containing the simulated points; a list of datasets (unless drop=TRUE and nsim=1).  Each dataset is a matrix with 2 or 3 columns, or an object of class sp2 or sp3
## Arguments:
## n : number of points to be simulated (numeric of length 1, > 0)
## mode : location parameter, gives the location of a mode (Watson bipolar) or antimode (Watson girdle).  Used to determine which simulation function to call.
## kappa : concentration parameter, numeric of length 1, must be > 0
## win : the window in which to simulate the points
## squeeze : logical.  If TRUE, then uses a squeeze function to accelerate the simulation process.  Recommended only when kappa is large.
## nsim : number of datasets to simulate.  numeric of length 1
## drop : logical.  If TRUE and nsim=1, the output will not be a list.  Ignored if nsim !=1
## as.sp : logical.  If TRUE, the output is a (list of) sp2 or sp3 objects, depending on the value of ndim
## ndim : numeric, equal to 2 or 3.  This determines the locations of points should be given in 2 dimensions (i.e. polar coordinates) or 3 dimensions (i.e. Cartesian coordinates)

rWatson.bipolar <- function(n, mode, kappa, win=sphwin(type="sphere"), squeeze=FALSE, nsim=1, drop=TRUE, as.sp=TRUE, ndim="2") {
	## initial checks of arguments
	stopifnot(inherits(win, "sphwin"))
	stopifnot(kappa > 0)
	stopifnot(sum(ndim==c("2", "3"))==1)
	stopifnot(n > 0)
	output <- list()
	## we now simulate each set of n points in turn
	for(i in 1:nsim) {
		X <- matrix(ncol=2, nrow=0)
		c1 <- 1/(exp(kappa)-1)
		## Having set up preliminaries, we now generate points 1 at a time
		while(nrow(X) < n) {
			u <- runif(2,0,1)
			u1 <- u[1]
			u2 <- u[2]
			s1 <- (1/kappa)*log((u1/c1)+1)
			test0 <- kappa*s1*(s1-1)
			## The first if argument applies the squeeze test if required.  If the squeeze test fails or is not required, the second if argument performs the standard test.  If either test passes, the point is will be added to dataset i (if addWatson confirms that it is in the window win)
			if(squeeze && (test0 + 1 >= u2)) {
				X <- addpoint.Watson(sround(cround(s1)), X, win, mode)
			} else if(exp(kappa*s1*(s1-1)) >= u2) {
 				X <- addpoint.Watson(sround(cround(s1)), X, win, mode)
			}
		}
		## Dataset i contains n points from a standard Watson bipolar distribution (i.e. mode=c(0,0)).  We now rotate this to the desired mode and take the required action for the output to be as defined by nsim, as.sp, ndim and (if applicable) drop
		X <- switch(ndim,
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
				output[[i]] <- sp2(X, win)
			}
			else {
				output[[i]] <- sp3(X, win)
			}
		}
	}
	if(nsim==1 && drop==TRUE) {
		output <- output[[1]]
	}
	output
}

## Function: rWatson.girdle
## Role: simulate datasets of n points from a Watson girdle distribution with location parameter mode and concentration parameter kappa. 
## Output: An object containing the simulated points; a list of datasets (unless drop=TRUE and nsim=1).  Each dataset is a matrix with 2 or 3 columns, or an object of class sp2 or sp3
## Arguments:
## n : number of points to be simulated (numeric of length 1, > 0)
## mode : location parameter, gives the location of a mode (Watson bipolar) or antimode (Watson girdle).  Used to determine which simulation function to call.
## kappa : concentration parameter, numeric of length 1, must be < 0
## win : the window in which to simulate the points
## squeeze : logical.  If TRUE, then uses a squeeze function to accelerate the simulation process.  Recommended only when kappa is large.
## nsim : number of datasets to simulate.  numeric of length 1
## drop : logical.  If TRUE and nsim=1, the output will not be a list.  Ignored if nsim !=1
## as.sp : logical.  If TRUE, the output is a (list of) sp2 or sp3 objects, depending on the value of ndim
## ndim : numeric, equal to 2 or 3.  This determines the locations of points should be given in 2 dimensions (i.e. polar coordinates) or 3 dimensions (i.e. Cartesian coordinates)


rWatson.girdle <- function(n, mode, kappa, win=sphwin(type="sphere"), squeeze=FALSE, nsim=1, drop=TRUE, as.sp=TRUE, ndim="2") {
	## initial checks of arguments
	stopifnot(inherits(win, "sphwin"))
	stopifnot(kappa < 0)
	stopifnot(sum(ndim==c("2", "3"))==1)
	stopifnot(n > 0)
	output <- list()
	## we now simulate each set of n points in turn
	for(i in 1:nsim) {
		X <- matrix(ncol=2, nrow=0)
		c1 <- sqrt(-kappa)
		c11 <- 1/c1
		c2 <- atan(c1)
		## Having set up preliminaries, we now generate points 1 at a time
		while(nrow(X) < n) {
			u <- runif(2,0,1)
			u1 <- u[1]
			u2 <- u[2]
			s1 <- c11*tan(c2*u1)
			test1 <- kappa*(s1^2)
			## The first if argument applies the squeeze test if required.  If the squeeze test fails or is not required, the second if argument performs the standard test.  If either test passes, the point is will be added to dataset i (if addWatson confirms that it is in the window win)
			if(squeeze && (1-(test1^2)) >= u2) {
				X <- addpoint.Watson(sround(cround(s1)), X, win, mode)
			} else if((1-test1)*exp(test1) >= u2) {
 				X <- addpoint.Watson(sround(cround(s1)), X, win, mode)
			}
		}
		## Dataset i contains n points from a standard Watson bipolar distribution (i.e. mode=c(0,0)).  We now rotate this to the desired mode and take the required action for the output to be as defined by nsim, as.sp, ndim and (if applicable) drop
		X <- switch(ndim,
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
				output[[i]] <- sp2(X, win)
			}
			else {
				output[[i]] <- sp3(X, win)
			}
		}
	}
	if(nsim==1 && drop==TRUE) {
		output <- output[[1]]
	}
	output
}

## Function: addpoint.Watson
## Internal function to add an observation from a standard Watson distribution to a set of such observations
## Arguments: 
## s1 : is the information required to determine the colatitude (theta) of the observation (the longitude phi is generated from a uniform distribution)
## X  : is the dataset to which the new point is added (a matrix with 2 columns)
## mode : the coordinates of a mode of the Watson distribution (matrix with 1 row and 2 or 3 columns, or numeric of length 2 or 3)
## win : is the window in which the points are to be located (the point defined by s1 is not added to X if it is outside win


addpoint.Watson <- function(s1, X, win, mode) {
	stopifnot(length(s1)==1 && s1 >= -1 && s1 <= 1 && inherits(X, "matrix") && ncol(X)==2 && inherits(win, "sphwin"))
	## Set the initial values of theta and phi (places the point in the upper hemisphere)
	twopi <- 2*pi
	theta <- acos(s1)
	phi <- runif(1,0,twopi)
	## Now perform a test to determine if the point should be in the lower hemisphere and if so, alter theta accordingly
	u3 <- runif(1,0,1)
	if(u3 >= 0.5) {
		theta <- pi-theta
	}
	point <- rot.sphere(c(theta, phi), northpole=mode)
	## Add the point, if it's in W
	if(in.W(matrix(point, nrow=1, ncol=2, byrow=TRUE), win)) {
		X <- rbind(X, point)
	}
	X
}