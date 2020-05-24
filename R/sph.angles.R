sph.angles <- function (win) {
	stopifnot(inherits(win, c("sp2", "sp3", "sphwin")))
	if(inherits(win, "sphwin")) {
		win <- win
	}
	else {
		win <- win$win
	}
	angles <- switch(win$type,
		wedge = {
			rep(win$param[1], 2)
		},
		quadrangle = {
			rep(pi/2, 4)
		},
		polygon = {
			sph.angles.poly(win=win)
		},
		{
		"This shape does not have any angles."
		}	
	)
	angles
}

## sph.angles.poly gives the algorithm for calculating the size of the angles in a polygon
## The algorithm being used here is split the angle at each vertex into two, by drawing a great circle arc from ref3 to the vertex that is entirely within the polygon.
## We can calculate the size of the two angles created to get the size of the angle at that vertex of the polygon (the sphcos commands at the end of the code; the preceding lines of code are used to get the information that the sphcos commands require)
## The algorithm is designed to simultaneously perform this calculation for all polygons.

sph.angles.poly <- function(win) {
	stopifnot(inherits(win, c("sp2", "sp3", "sphwin")))
	if(inherits(win, "sphwin")) {
		win <- win
	}
	else {
		win <- win$win
	}
	stopifnot(win$type=="polygon")
	rad <- win$rad
	param3 <- convert3(win$param)
	lp <- nrow(param3)
	ref3 <- matrix(convert3(win$ref3), nrow=1, ncol=3, byrow=TRUE)
	polyangs1 <- param3[1:(lp-1),]
	polyangs2 <- param3[2:lp,]
	polyangs3 <- rbind(param3[3:lp, ], param3[2,])
	gc <- gcdistPaired(x=polyangs1, y=polyangs2, rad=rad)
	gc1 <- gc[1:(lp-1)]
	gc2 <- c(gc[2:(lp-1)], gc[1])
	gc.r3v1 <- gcdist(x=ref3, y=polyangs1)
	gc.r3v2 <- gcdist(x=ref3, y=polyangs2)
	gc.r3v3 <- gcdist(x=ref3, y=polyangs3)
	sc1 <- sphcos(d1=gc1, d2=gc.r3v2, d3=gc.r3v1, theta=NULL, rad=rad)
	sc2 <- sphcos(d1=gc2, d2=gc.r3v2, d3=gc.r3v3, theta=NULL, rad=rad)
	angles <- sc1+sc2	
	angles
}
