## is.sphcoords - tests whether a point x is on the sphere with radius rad, and is in spherical coordinates
## x needs to be a matrix with one row and two columns or a numeric of length 2 in order to pass the initial test
## rad is a non-negative length 1 numeric
## Internal function used in sphwin

is.sphcoords <- function(x, rad) {	
	stopifnot(rad > 0)
	issphcoords <- FALSE
	if(inherits(x, "matrix")) {
		if(nrow(x)==1 && ncol(x)==2) {
			x1 <- x[1,1]/rad
			x2 <- x[1,2]/rad	
			if(x1 >= 0 && x1 <= pi && x2 >=0 && x2 < 2*pi) {
				issphcoords <- TRUE
			}
		}
	} else if(length(x)==2) {
		x1 <- x[1]/rad
		x2 <- x[2]/rad	
		if(x1 >= 0 && x1 <= pi && x2 >=0 && x2 < 2*pi) {
			issphcoords <- TRUE
		}
	}
	issphcoords	
}