## Function: in.W
## Determines whether a point is (or points are) in a window
## Arguments: points - a 2 or 3 column matrix containing the locations of points
##	      win - the window, an object of class sphwin
## Output: a length n logical vector (n being the number of points)
##         with ith entry TRUE if point i is in W,
##         and FALSE if point i is not in W.


in.W <- function(points, win) {

  ## Preliminary check to ensure win is an object of class sphwin
  stopifnot(inherits(win, "sphwin"))

  X <- if(is.matrix(points)) points else
       if(inherits(points, c("sp2", "sp3"))) points$X else
       stop("points should be a matrix or an object of class sp2 or sp3")

  nX <- nrow(X)
  if(nX == 0) return(logical(0))

  rad <- win$rad
  radpar <- rad * win$param
  
  switch(win$type,
         sphere = {
           return(rep(TRUE, nX))
         },
         band = {
           centre <- matrix(win$ref, nrow=1)
           dcentre <- as.vector(gcdist(X, centre, rad=rad))
           result <- (radpar[1] <= dcentre) & (dcentre <= radpar[2])
           return(result)
         },
         bandcomp = {
           centre <- matrix(win$ref, nrow=1)
           dcentre <- as.vector(gcdist(X, centre, rad=rad))
           result <- (dcentre <= radpar[1]) | (dcentre >= radpar[2])
           return(result)
         },
         polygon = {
           result <- in.W.poly(points=points, win=win)
           return(result)
         })

  ## The window is a wedge or quadrangle.
  ## We first rotate the points so that the window's boundaries
  ## are lined up with circles of colatitude (quadrangle)
  ## and semicircles of longitude (wedge, quadrangle).
  ## This simplifies significantly the check we need to do
  ## to determine whether the point or points are in the window.

  if(identical(win$ref, c(0,0))) {
    points1 <- points
  } else {
    points1 <- rot.sphere(points=points, northpole=win$ref,
                          rad=rad, inverse=FALSE)
  }

  ## If the window is a wedge, then the previous rotation may result
  ## in the left-hand edge of the wedge not being at longitude 0,
  ## so we rotate again to achieve that characteristic.
  ## This simplifies significantly the check we need to do
  ## to determine whether the point or points are in the window.

  if(win$type=="wedge") {
    points1 <- rot.sphere(points=points1, northpole=c(0,win$param[2]),
                          rad=rad, inverse=FALSE)
  }
	
  ## Next, we perform the relevant check:
  ## for the quadrangle we need to check whether the
  ## colatitudes of points are within the correct range, ditto for
  ## longitudes of the rotated points for the wedge and quadrangle. 

  #' hack inserted by Adrian
  if(inherits(points1, c("sp2", "sp3")))
    points1 <- points1$X
  #' end hack
           
  m <- nrow(points1)
  colats <- points1[,1]
  lons <- points1[,2]
  winparam1 <- win$param[1]
  winparam2 <- win$param[2]
  is.in.W <-
    switch(win$type,
           wedge={
             (sround(lons)>=0) & (sround(winparam1-lons)>=0)
           },
           quadrangle = {
             winparam3 <- win$param[3]
             inWq <- (sround(colats - winparam1) >= 0) &
                     (sround(winparam2 - colats) >= 0) &
                     (sround(lons) >= 0) &
                     (sround(winparam3-lons) >=0)
             inWq
           },
           {stop("Unrecognised window type")}
           )
  is.in.W
}

