## Function: in.W
## Determines whether a point is (or points are) in a window
## Arguments: points - a 2 or 3 column matrix containing the locations of points
##	      win - the window, an object of class sphwin
## Output: a length n numeric (n being the number of points) with ith entry TRUE if point i is in W, and FALSE if point i is not in W.


in.W <- function(points, win) {

	## Preliminary check to ensure win is an object of class sphwin

	stopifnot(inherits(win, "sphwin"))

	## If the window is a band complement, a point in the window is in one of the two caps that form the window, so we make a fresh call to in.W accordingly:

	if(win$type=="bandcomp") {
		points <- rot.sphere(points=points, northpole=win$ref, inverse=FALSE)
		inWs <- in.W(points, sphwin(type="band", param=c(0, win$param[1]), ref=c(0,0))) + in.W(points, sphwin(type="band", param=c(win$param[2], pi), ref=c(0,0)))
		is.in.W <- sum(inWs >= 1)
	} else {
	
	## If the window is a band, wedge or quadrangle, we first rotate the points so that the window's boundaries are lined up with circles of colatitude (band, quadrangle) and semicircles of ## longitude (wedge, quadrangle).  This simplifies significantly the check we need to do to determine whether the point or points are in the window.

	if(win$type=="polygon" || win$type=="sphere" || identical(win$ref, c(0,0))) {points1 <- points} else { points1 <- rot.sphere(points=points, northpole=win$ref, inverse=FALSE)}
	rad <- win$rad

	## If the window is a wedge, then the previous rotation may result in the left-hand edge of the wedge not being at longitude 0, so we rotate again to achieve that characteristic.  This ## simplifies significantly the check we need to do to determine whether the point or points are in the window.

	if(win$type=="wedge") {points1 <- rot.sphere(points=points1, northpole=c(0,win$param[2]), inverse=FALSE)}
	
	## Next, we perform the relevant check: for the sphere no checks are needed, for the band and quadrangle we need to check whether the colatitudes of points are within the correct range, ditto for longitudes of the rotated points for the wedge and quadrangle.  The polygon is a more complicated situation - the function in.W.poly performs the check for us.

        #' hack inserted by Adrian
        if(inherits(points1, c("sp2", "sp3")))
          points1 <- points1$X
        #' end hack
           
	m <- nrow(points1)
	colats <- points1[,1]
	lons <- points1[,2]
	winparam1 <- win$param[1]
	winparam2 <- win$param[2]
	is.in.W <- switch(win$type,
	sphere= {rep(TRUE,m)},
	band={(sround(colats - winparam1) >= 0) + (sround(winparam2 - colats) >= 0) ==2},
	wedge={(sround(lons)>=0) + (sround(winparam1-lons)>=0)==2},
	polygon ={in.W.poly(points=points, win=win)},
	quadrangle = {
		      winparam3 <- win$param[3]
		      inWq <- ((sround(colats - winparam1) >= 0) + (sround(winparam2 - colats) >= 0) + (sround(lons)>=0) + (sround(winparam3-lons)>=0) ==4)
		      inWq
		     },
{stop("Unrecognised window type")}
)

}
is.in.W
}

