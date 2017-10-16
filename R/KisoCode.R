## File KisoCode.R
##   (was: '20140722 Kiso code.txt')

## Function: intcircs
## intcircs calculates the intersection points between two circles on a sphere.  This covers everything in Chapter 0.1.1
## n1 is a vector normal to the first circle.  It must be a 3D vector.
## h1 is the dot product of the vector between the centre of the sphere and the first circle, and vector between the centre of the sphere and any point on the circumference of the  first circle
## h1 must be length 1.
## n2 and h2 are the corresponding values for the second circle.
## n1 and h1 are such that the first circle can be defined as {x \in r^3 : x.n1 = h1} where x is a point on the sphere and . is the dot product.
## n2 and h2 take the corresponding meaning for the second circle.
## tol is the tolerance for deciding when lambda2==0 (and hence lambda==0)

.InfiniteSet <- matrix(Inf, nrow=1, ncol=1)
.EmptySet    <- matrix(NA, nrow=1, ncol=1)

intcircs <- function(n1, n2, h1, h2, tol=1e-12) {
  ## Validate the input data
  stopifnot(length(n1) == 3)
  stopifnot(length(n2) == 3)
  stopifnot(length(h1)==1 && length(h2)==1)
	
  ## Form cross product
  cn <- sround(cross(n1, n2))

  ## The squared modulus of cross product
  scn2 <- sround(sum(cn^2))

  ## The dot product of the normals
  dn <- sround(dot(n1, n2))

  ## Now the cross product, its modulus and the dot product
  ## have been calculated, we find c1, c2 (see equation above (0.3) )
  ## and from that we can find lambda2 (=lambda^2)
  ## We first consider only cases where the length of the cross product
  ## is nonzero

  if(scn2==0) {
    ## Intersection is either empty or infinite;
    ## output is a 1 cell matrix, containing Inf if the intersection
    ## is infinite or NA if the intersection is empty.
    ## Either n1==n2 or n1==-n2
    if((dn > 0 && abs(h1-h2) < 1e-12) || (dn < 0 && abs(h1+h2) < 1e-12) ) {
      output <- .InfiniteSet
    } else {
      output <- .EmptySet
    }
  } else {
    ## Finite nonzero number of intersection points
    c1 <- sround((h1 - (h2*dn))/scn2)
    c2 <- sround((h2 - (h1*dn))/scn2)
    const <- sround((c1*n1) + (c2*n2))

    ## The next line is equation (0.7) in the appendix
		
    t2 <- sround((1-(c1^2)-(c2^2)-(2*c1*c2*dn))/scn2)
		
    ## To avoid numerical problems, we treat lambda <= 10^-4
    ## as being lambda ==0.  If lambda2 >= 0 we can find the
    ## intersection point(s)

    if(t2 >= tol) {
      t1 <- sqrt(t2)*c(-1, 1)
      output <- rbind(sround(const + t1[1]*cn), sround(const+t1[2]*cn))
    } else if(t2 < tol) {
      output <- matrix(const, nrow=1, ncol=3, byrow=TRUE)
    } else 	
      output <- .EmptySet
  }
  return(output)
}

# Functions: Kiso
# Calculates Kiso for the wedge, band/cap and polygon respectively (called by Kiso, rather than directly by the user)
# Arguments: X is an object of class sp2 or sp3, or a two or three column matrix giving the spherical or cartesian coordinates of points in the point pattern
#            win is an object of type sphwin
#            r is the vector of distances at which the estimator is calculated


Kiso <- function(X, win, r, rad=win$rad, Dmat=pairdistsph(X), nrX=nrow(X), denom=(nrX*(nrX-1))/area.sphwin(X), lambda=NULL) {

  ## First, we ensure the variables are of the acceptable input type,
  ## and if X is of class sp2 or sp3, we save the locations of points,
  ## and the window as separate objects 

  stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
  if(inherits(X, "matrix")) {
    stopifnot(inherits(win, "sphwin"))
  } else {
    win <- X$win
    X <- X$X
  }

  ## For efficiency, we save the radius as a separate object,
  ## and ensure the locations of points are in spherical coordinates

  rad <- win$rad
  if(ncol(X)==3) {
    X <- convert2(points=X, rad=rad)
  }
  type <- win$type

  ## We now calculate wmat (i.e. equation (0.1)
  wmat <- switch(type,
                 sphere = matrix(1, nrX, nrX),
                 band = {
                   if(win$param[1]==0 || win$param[2]==pi) {
                     Kisocap(X=X, win=win, r=r, nrX=nrX, Dmat=Dmat, disc=FALSE,
                             denom=denom, lambda=lambda)
                   } else {
                     Kisoband(X=X, win=win, r=r, nrX=nrX, Dmat=Dmat, disc=FALSE,
                              denom=denom, lambda=lambda)
                   }
                 },
                 bandcomp = {
                   Kisobc(X=X, win=win, r=r, nrX=nrX, Dmat=Dmat, disc=FALSE,
                          denom=denom, lambda=lambda)
                 },
                 wedge = {
                   Kisowedge(X=X, win=win, r=r, nrX=nrX, Dmat=Dmat, disc=FALSE)
                 },
                 polygon = {
                   Kisopoly(X=X, win=win, r=r, nrX=nrX, Dmat=Dmat, disc=FALSE)
                 },
		 quadrangle = {
                   Kisoquad(X=X, win=win, r=r, nrX=nrX, Dmat=Dmat, disc=FALSE)
                 },
                 stop("Unrecognised shape type")
                 )

  if(!is.null(lambda) &&
     type %in% c("sphere", "wedge", "polygon", "quadrangle")) 
    wmat <- wmat * outer(lambda, lambda, "*")
    
  ## shortcut
  if(is.fv(wmat)) return(wmat)

  ## otherwise wmat is a matrix
  if(any(is.na(wmat))) {stop("NA in weight matrix")}

  ## Finally, we can make a call to compileK to calculate Kiso
  Kisoout <- compileK(D=Dmat, r=r, weights=1/wmat, denom=denom)
  return(Kisoout)
}


## Function: Kisocap
## Calculate: Kiso for the cap (called by Kiso, rather than directly by the user)
## Arguments: X is a two column matrix giving the spherical coordinates of points in the point pattern
##            win is an object of type sphwin
##            r is the vector of distances at which the estimator is calculated
##            nrX, Dmat are respectively the number of rows in X, and the pairwise distance matrix for X
##            rad is the radius of the sphere (it is an element of sphwin but having it separately below is more efficient)
##            disc is a logical, if it is TRUE then the discretizes estimation of the weight matrix is performes]d

Kisocap <- function(X, win, r, nrX=nrow(X), Dmat=pairdistsph(X),
                    disc=FALSE, rad=win$rad, denom, useC, lambda=NULL) {

  if(missing(useC)) useC <- getOption("sphwin.useC")
  if(is.null(useC)) useC <- TRUE

  ## Preliminary checks, and setting up objects that are required later
  stopifnot(win$type=="band" && (win$param[1]==0 || win$param[2]==pi))
  if(win$param[1]==0) {
    caprad <- win$param[2]
  } else {
    caprad <- win$param[1]
  }
  capheight <- cos(caprad)

  ## new C code
  if(useC) {
    centre <- convert3(win$ref)
    if(ncol(X) != 3) X <- convert3(X)
    n <- nrow(X)
    if(!is.null(Dmat) && is.null(lambda)) {
      ## calculate weight matrix
      zz <- .C("kisocapweights",
               n = as.integer(n),
               x1 = as.double(X[,1]),
               x2 = as.double(X[,2]),
               x3 = as.double(X[,3]),
               Dmat = as.double(Dmat),
               centre = as.double(centre),
               height = as.double(capheight),
               wmat = as.double(numeric(n * n)))
      return(matrix(zz$wmat, n, n))
    } else if(is.null(lambda)) {
      ## calculate empirical K-function
      nr <- length(r)
      rmax <- max(r)
      zz <- .C("dkisocap",
               n = as.integer(n),
               x1 = as.double(X[,1]),
               x2 = as.double(X[,2]),
               x3 = as.double(X[,3]),
               centre = as.double(centre),
               height = as.double(capheight),
               nr = as.integer(nr),
               rmax = as.double(rmax),
               dk = as.double(numeric(nr)))
      kval <- cumsum(zz$dk)/denom
      df <- data.frame(r=r, est=kval)
      K <- fv(df, "r", quote(K(r)), "est", . ~ r, c(0, rmax), 
            c("r", "hat(%s)[iso](r)"),
              c("distance argument r", 
                "isotropic correction estimate of %s"),
              fname = "K")
      return(K)
    } else {
      ## calculate inhomogeneous K-function
      check.nvector(lambda, n)
      nr <- length(r)
      rmax <- max(r)
      zz <- .C("dwkisocap",
               n = as.integer(n),
               x1 = as.double(X[,1]),
               x2 = as.double(X[,2]),
               x3 = as.double(X[,3]),
               lambda = as.double(lambda),
               centre = as.double(centre),
               height = as.double(capheight),
               nr = as.integer(nr),
               rmax = as.double(rmax),
               dk = as.double(numeric(nr)))
      kval <- cumsum(zz$dk)/denom
      df <- data.frame(r=r, est=kval)
      K <- fv(df, "r", quote(K(r)), "est", . ~ r, c(0, rmax), 
            c("r", "hat(%s)[iso](r)"),
              c("distance argument r", 
                "isotropic correction estimate of %s"),
              fname = "K")
      return(K)
    }
  }
  
  ## Initialise all weights to 1
  wmat <- matrix(1, nrow=nrX, ncol=nrX)
  diag(wmat) <- 1
	
  Drad <- Dmat/rad
  CD <- sround(cround(cos(Drad)))
  winref <- win$ref
  winref3 <- convert3(win$ref)

  ## X is in spherical coordinates; convert to Cartesian coordinates
  X3 <- convert3(X)

  if(disc) {
    ## Also calculate the discretized estimator.
    circs.disc <- matrix(nrow=nrX, ncol=nrX)
    diag(circs.disc) <- 1
    lons <- seq(0, 2*pi*(1-(100^-1)), length=100)

    for(i in 1:nrX) {
      ## Create objects for X[i,] in spherical and Cartesian coordinates
      xi <- X[i,]
      xi3 <- X3[i,]
      for(j in setdiff(1:nrX, i)) {
        ## Create objects for X[j,] in spherical and Cartesian coordinates
        xj <- X[j,]
        xj3 <- X3[j, ]
        ## cosine of distance between xi and xj
        cdij <- CD[i,j]
        ## Calculate the intersection points between the cap and \partial bij
        ints <- intcircs(xi3, winref3, cdij, capheight)
      
        ## ..... redundant ..........
        ## Cases C1, C2 and C3
        ## if(identical(ints, .InfiniteSet) ||
        ##   identical(ints, .EmptySet) ||
        ##   nrow(ints)==1) {
        ##  wmat[i,j] <- 1
        ## } else {
        ## ............................
      
        ## Case C4
        if(nrow(ints)==2) {
          wmat[i,j] <- Kisoengine(xi3=xi3, xj3=xj3, win=win,
                                  ints=ints, cdij=cdij)
        }

        ints.disc.ij <- rot.sphere(cbind(rep(Drad[i,j], times=100), lons),
                                   northpole=xi, inverse=TRUE)
        circs.disc[i,j] <- sum(in.W(points=ints.disc.ij, win=win))/100
      }
    }
    ## The code below turns the output into a list containing the
    ## weights matrices (each of which could be provided to compileK
    ## the argument \code{weights}) for the actual and discretized estimators
    ## of Kiso
    if(disc) {
      attr(wmat, "discrete") <- circs.disc
    }
  } else {
    ## SAME AS ABOVE, BUT WITHOUT DISCRETISED ESTIMATES
    for(i in 1:nrX) {
      ## Create objects for X[i,] in spherical and Cartesian coordinates
      xi <- X[i,]
      xi3 <- X3[i,]
      for(j in setdiff(1:nrX, i)) {
        ## Create objects for X[j,] in spherical and Cartesian coordinates
        xj <- X[j,]
        xj3 <- X3[j, ]
        ## cosine of distance between xi and xj
        cdij <- CD[i,j]
        ## Calculate the intersection points between the cap and \partial bij
        ints <- intcircs(xi3, winref3, cdij, capheight)
        ## Case C4
        if(nrow(ints)==2) {
          wmat[i,j] <- Kisoengine(xi3=xi3, xj3=xj3, win=win,
                                  ints=ints, cdij=cdij)
        }
      }
    }
  }
  if(anyNA(wmat)) {stop("NA in weight matrix")}
  return(wmat)
}

## Function: Kisoengine
## Calculate: wij/||b_ij|| given the intersection points, xi and xj
## Arguments: xi3 is the centre of bij, given in Cartesian coordinates
##	      xj3 is the point in the point pattern that is on the boundary of bij, in Cartesian coordinates
##            win is an object of type sphwin
##	      ints is a matrix containing the intersection points, in Cartesian coordinates

Kisoengine <- function(xi3, xj3, win, ints, verbose=FALSE, cdij=dot(xi3, xj3)) {
  ## define the orthonormalised vectors wi, vi, si as in (0.5)
  wi  <- xj3 - cdij*xi3
  vi <- wi/sqrt(sum(wi^2))
  si <- cross(xi3, vi)
#  nri <- nrow(ints) # never used

  if(verbose) {
    cat("ints=\n"); print(ints);
    cat("vi=\n"); print(vi);
    cat("si=\n"); print(si)
  }

  ## find the intersection points in the coordinates of the orthonormalization,
  ## and in the final angular coordinate

  thetaint <- c()

  ## This matrix below contains the u.vi and u.si values defined in (0.7),
  ## and theta as defined in (0.8)

  visi <- cbind(ints %*% matrix(vi, ncol=1, nrow=3),
                ints %*% matrix(si, ncol=1, nrow=3))

  if(verbose) { cat("visi=\n"); print(visi) }
  thetaint <- atan2(visi[,2], visi[,1]) %% (2*pi)
  nints <- length(thetaint)

  if(verbose) { cat("thetaint=\n"); print(thetaint) }

  ## If nints=1, we need a second intersection point,
  ## so we include the point that is diametrically opposite
  ## the point thetaint in the cap

  ## We order the intersection points, and put the first one
  ## at the end of the vector for use in the code that follows

  sortints <- sort(thetaint)
  sortints <- c(sortints, sortints[1] + 2*pi)
  
  if(verbose) { cat("sortints=\n"); print(sortints) }
	
  ## We find the points halfway between each pair of intersection points
  ## in terms of the angular coordinate, then the orthonormalized coordinates

  if(nints > 1) {
    alphas <- ((sortints[2:(nints+1)] + sortints[1:nints])/2) %% (2*pi)
    sortintdiffs <- diff(sortints)
  } else if (nints==1) {
    alphas <- (thetaint + pi) %% (2*pi)
    sortintdiffs <- 1
  }
  xdij <- xi3*as.numeric(cdij)

  if(verbose) { cat("alphas=\n"); print(alphas) }
  if(verbose) { cat("sortintdiffs=\n"); print(sortintdiffs) }

  wij <- 0
  sdij <- sqrt(1 - cdij^2)
  for(j in seq_len(nints)) {
    ## find the Cartesian coordinate of each midpoint,
    ## test whether it is in W, and if a midpoint is in W
    ## add the length of its arc to wij.  This corresponds to equation (0.9).

#    qadd <-  sqrt(1-(cdij^2))*(vi*cos(alphas[j]) + si*sin(alphas[j])) + xdij
    qadd <-  sdij * (vi*cos(alphas[j]) + si*sin(alphas[j])) + xdij
		
    if(verbose) { cat("qadd=\n"); print(qadd) }

    if(in.W(convert2(qadd), win)) {
      
      if(verbose) { cat("in.W="); print(TRUE) }

      wij <- wij+sortintdiffs[j]
    } else { if(verbose) { cat("in.W="); print(FALSE) } }
  }

  ## Normalize wij (i.e. convert it to wij)	
  wij <- wij/(2*pi)
  return(wij)
}


## Function: Kisoband
## Calculate: Kiso for the band (called by Kiso, rather than directly by the user)
## Arguments: X is a two column matrix giving the spherical coordinates of points in the point pattern
##            win is an object of type sphwin
##            r is the vector of distances at which the estimator is calculated
##            nrX, Dmat are respectively the number of rows in X, and the pairwise distance matrix for X
##            rad is the radius of the sphere (it is an element of sphwin but having it separately below is more efficient)
##            disc is a logical, if it is TRUE then the discretized estimation of the weight matrix is performed



Kisoband <- function(X, win, r, nrX=nrow(X), Dmat=pairdistsph(X), disc=FALSE, rad=win$rad, denom, useC, lambda=NULL) {

  if(missing(useC)) useC <- getOption("sphwin.useC")
  if(is.null(useC)) useC <- TRUE
        
  ## Checking that win is actually a band,
  ## and making the bounding circles of colatitude objects

  stopifnot(win$type=="band" && win$param[1]!=0 && win$param[2]!=pi)
  lat1 <- win$param[1]
  lat2 <- win$param[2]
  clat1 <- cos(lat1)
  clat2 <- cos(lat2)
  winref <- win$ref

  if(useC) {
    centre <- convert3(winref)
    if(ncol(X) != 3) X <- convert3(X)
    n <- nrow(X)
    if(!is.null(Dmat) && is.null(lambda)) {
      ## calculate weight matrix
      zz <- .C("kisobandweights",
               n = as.integer(n),
               x1 = as.double(X[,1]),
               x2 = as.double(X[,2]),
               x3 = as.double(X[,3]),
               Dmat = as.double(Dmat),
               centre = as.double(centre),
               height1 = as.double(clat1),
               height2 = as.double(clat2),
               iscomp = as.integer(0),
               wmat = as.double(numeric(n * n)))
      return(matrix(zz$wmat, n, n))
    } else if (is.null(lambda)) {
      ## calculate empirical K-function
      nr <- length(r)
      rmax <- max(r)
      zz <- .C("dkisoband",
               n = as.integer(n),
               x1 = as.double(X[,1]),
               x2 = as.double(X[,2]),
               x3 = as.double(X[,3]),
               centre = as.double(centre),
               height1 = as.double(clat1),
               height2 = as.double(clat2),
               iscomp = as.integer(0),
               nr = as.integer(nr),
               rmax = as.double(rmax),
               dk = as.double(numeric(nr)))
      kval <- cumsum(zz$dk)/denom
      df <- data.frame(r=r, est=kval)
      K <- fv(df, "r", quote(K(r)), "est", . ~ r, c(0, rmax), 
            c("r", "hat(%s)[iso](r)"),
              c("distance argument r", 
                "isotropic correction estimate of %s"),
              fname = "K")
      return(K)
    } else {
      ## calculate inhomogeneous K-function
      check.nvector(lambda, n)
      nr <- length(r)
      rmax <- max(r)
      zz <- .C("dwkisoband",
               n = as.integer(n),
               x1 = as.double(X[,1]),
               x2 = as.double(X[,2]),
               x3 = as.double(X[,3]),
               lambda = as.double(lambda),
               centre = as.double(centre),
               height1 = as.double(clat1),
               height2 = as.double(clat2),
               iscomp = as.integer(0),
               nr = as.integer(nr),
               rmax = as.double(rmax),
               dk = as.double(numeric(nr)))
      kval <- cumsum(zz$dk)/denom
      df <- data.frame(r=r, est=kval)
      K <- fv(df, "r", quote(K(r)), "est", . ~ r, c(0, rmax), 
            c("r", "hat(%s)[iso](r)"),
              c("distance argument r", 
                "isotropic correction estimate of %s"),
              fname = "K")
      return(K)
    }
  }
        
  ## Extract the two bounding circles as separate windows

	win1 <- sphwin(type="band", param=c(lat1, pi), ref=winref)
	win2 <- sphwin(type="band", param=c(0, lat2), ref=winref)
	
	## If we want the discretized estimator, this calls that object	

	if(disc) {
			circs.disc <- matrix(nrow=nrX, ncol=nrX)
			diag(circs.disc) <- 1
			lons <- seq(0, 2*pi*(1-(100^-1)), length=100)
	}

	## Creates the output object

	wmat <- matrix(nrow=nrX, ncol=nrX)
	diag(wmat) <- 1
	Drad <- Dmat/rad
	CD <- sround(cround(cos(Drad)))
	winref <- win$ref
	winref3 <- convert3(winref)
	for(i in 1:nrX) {

		## Get data for xi in both spherical and Cartesian coordinates

		xi <- X[i,]
		xi3 <- convert3(xi)
		for(j in 1:nrX) {
			if(i != j) {

				## Get data for xj in both spherical and Cartesian coordinates

				xj <- X[j,]
				xj3 <- convert3(xj)

				## Get data for the distance between and dot product of xi and xj

				dij <- Dmat[i,j]
# ajb: 'dij' is defined but never used                                
				cdij <- CD[i,j]

				## find all crossing points

				ints1 <- intcircs(n1=winref3, n2=xi3, h1=clat1, h2=cdij)
				ints2 <- intcircs(n1=winref3, n2=xi3, h1=clat2, h2=cdij)
				nr1 <- nrow(ints1)
				nr2 <- nrow(ints2)


				## If there are two crossing points with one of the circles of colatitude and either one or zero with the other, we calculate wij using only the circle of colatitude 				## that had two crossing points.  If there are two crossing points with both circles of colatitude then we calculate wij in a similar way.  Otherwise, wij=1.

				if(nr1 == 2 && nr2 < 2) {

					## Case B3 (2 crossing points in I1, less than 2 in I2)

					wmat[i,j] <- Kisoengine(xi3=xi3, xj3=xj3, win=win1, ints=ints1, cdij=cdij)				

				} else if(nr1 < 2 && nr2 == 2) {

					## Case B3 (2 crossing points in I2, less than 2 in I1)
					
					wmat[i,j] <- Kisoengine(xi3=xi3, xj3=xj3, win=win2, ints=ints2, cdij=cdij)
					
				} else if(nr1==2 && nr2==2) {

					## Case B4

					wmat[i,j] <- Kisoengine(xi3=xi3, xj3=xj3, win=win, ints=rbind(ints1, ints2), cdij=cdij)
					
				} else if(nr1 + nr2 <= 2) {
				
					## Cases B1 and B2

					wmat[i,j] <- 1
				}

				if(disc) {

					## Also calculate discretized approximation of the same quantity

					ints.disc.ij <- rot.sphere(cbind(rep(Drad[i,j], times=100), lons), northpole=xi, inverse=TRUE)
					circs.disc[i,j] <- sum(in.W(points=ints.disc.ij, win=win))/100
				}
			}
		}
	}

## In the event we want the discretized estimator calculated, the code below turns the output into a data.frame containing the weights matrices
## (each of which could be provided to compileK the argument \code{weights}) for the actual and discretizedt estimators of Kiso

	if(disc) {
		attr(wmat, "discrete") <- circs.disc
	}
	if(any(is.na(wmat))) {stop("NA in weight matrix")}
	return(wmat)
}

## Function: Kisobc
## Calculate: Kiso for the band (called by Kiso, rather than directly by the user)
## Arguments: X is a two column matrix giving the spherical coordinates of points in the point pattern
##            win is an object of type sphwin
##            r is the vector of distances at which the estimator is calculated
##            nrX, Dmat are respectively the number of rows in X, and the pairwise distance matrix for X
##            rad is the radius of the sphere (it is an element of sphwin but having it separately below is more efficient)
##            disc is a logical, if it is TRUE then the discretized estimation of the weight matrix is performed



Kisobc <- function(X, win, r, nrX=nrow(X), Dmat=pairdistsph(X), disc=FALSE, verbose=FALSE, rad=win$rad, denom, useC, lambda=NULL) {

  if(missing(useC)) useC <- getOption("sphwin.useC")
  if(is.null(useC)) useC <- TRUE
        
  ## Checking that win is actually a band complement,
  ## and making the bounding circles of colatitude objects

  stopifnot(win$type =="bandcomp" && win$param[1]!=0 && win$param[2]!=pi)
  lat1 <- win$param[1]
  lat2 <- win$param[2]
  clat1 <- cos(lat1)
  clat2 <- cos(lat2)
	
  winref <- win$ref
  winref3 <- convert3(winref)

  if(useC) {
    if(ncol(X) != 3) X <- convert3(X)
    n <- nrow(X)
    if(!is.null(Dmat) && is.null(lambda)) {
      ## compute matrix of correction weights
      zz <- .C("kisobandweights",
               n = as.integer(n),
               x1 = as.double(X[,1]),
               x2 = as.double(X[,2]),
               x3 = as.double(X[,3]),
               Dmat = as.double(Dmat),
               centre = as.double(winref3),
               height1 = as.double(clat1),
               height2 = as.double(clat2),
               iscomp = as.integer(1),
               wmat = as.double(numeric(n * n)))
      return(matrix(zz$wmat, n, n))
    } else if(is.null(lambda)) {
      ## compute empirical K-function
      nr <- length(r)
      rmax <- max(r)
      zz <- .C("dkisoband",
               n = as.integer(n),
               x1 = as.double(X[,1]),
               x2 = as.double(X[,2]),
               x3 = as.double(X[,3]),
               centre = as.double(winref3),
               height1 = as.double(clat1),
               height2 = as.double(clat2),
               iscomp = as.integer(1),
               nr = as.integer(nr),
               rmax = as.double(rmax),
               dk = as.double(numeric(nr)))
      kval <- cumsum(zz$dk)/denom
      df <- data.frame(r=r, est=kval)
      K <- fv(df, "r", quote(K(r)), "est", . ~ r, c(0, rmax), 
            c("r", "hat(%s)[iso](r)"),
              c("distance argument r", 
                "isotropic correction estimate of %s"),
              fname = "K")
      return(K)
    } else {
      ## compute inhomogeneous K-function
      check.nvector(lambda, n)
      nr <- length(r)
      rmax <- max(r)
      zz <- .C("dwkisoband",
               n = as.integer(n),
               x1 = as.double(X[,1]),
               x2 = as.double(X[,2]),
               x3 = as.double(X[,3]),
               lambda = as.double(lambda),
               centre = as.double(winref3),
               height1 = as.double(clat1),
               height2 = as.double(clat2),
               iscomp = as.integer(1),
               nr = as.integer(nr),
               rmax = as.double(rmax),
               dk = as.double(numeric(nr)))
      kval <- cumsum(zz$dk)/denom
      df <- data.frame(r=r, est=kval)
      K <- fv(df, "r", quote(K(r)), "est", . ~ r, c(0, rmax), 
            c("r", "hat(%s)[iso](r)"),
              c("distance argument r", 
                "isotropic correction estimate of %s"),
              fname = "K")
      return(K)
    }
  }
  

	## If we want the discretized estimator, this calls that object	

	if(disc) {
			circs.disc <- matrix(nrow=nrX, ncol=nrX)
			diag(circs.disc) <- 1
			lons <- seq(0, 2*pi*(1-(100^-1)), length=100)
			Drad <- Dmat/rad
	}

	## Creates the output object

	wmat <- matrix(nrow=nrX, ncol=nrX)
	diag(wmat) <- 1
	
	## Extract the two bounding circles as separate windows

	win1 <- sphwin(type="band", param=c(0, lat1), ref=winref)
	win2 <- sphwin(type="band", param=c(lat2, pi), ref=winref)
	Drad <- Dmat/rad
	CD <- sround(cround(cos(Drad)))
	for(i in 1:nrX) {

		## Get data for xi in both spherical and Cartesian coordinates

		xi <- X[i,]
		xi3 <- convert3(xi)
		for(j in 1:nrX) {
			if(i != j) {

				## Get data for xj in both spherical and Cartesian coordinates

				xj <- X[j,]
				xj3 <- convert3(xj)

				## Get data for the distance between and dot product of xi and xj

				dij <- Dmat[i,j]
# ajb: 'dij' is defined but never used                                
				cdij <- CD[i,j]

				## find all crossing points

				ints1 <- intcircs(n1=winref3, n2=xi3, h1=clat1, h2=cdij)
				ints2 <- intcircs(n1=winref3, n2=xi3, h1=clat2, h2=cdij)
				nr1 <- nrow(ints1)
				nr2 <- nrow(ints2)

				## If there are two crossing points with one of the circles of colatitude and either one or zero with the other, we calculate wij using only the circle of colatitude 				## that had two crossing points.  If there are two crossing points with both circles of colatitude then we calculate wij in a similar way.  Otherwise, wij=1.

				if(nr1 == 2 && nr2 < 2) {

					## Case BC4 (2 crossing points in I1, less than 2 in I2)

					wmat[i,j] <- Kisoengine(xi3=xi3, xj3=xj3, win=win1, ints=ints1, cdij=cdij, verbose)				

				} else if(nr1 < 2 && nr2 == 2) {

					## Case BC4 (2 crossing points in I2, less than 2 in I1)
					
					wmat[i,j] <- Kisoengine(xi3=xi3, xj3=xj3, win=win2, ints=ints2, cdij=cdij, verbose)
					
				} else if(nr1==2 && nr2==2) {

					## Case BC5 (2 crossing points in each of I1 and I2)

					wmat[i,j] <- Kisoengine(xi3=xi3, xj3=xj3, win=win, ints=rbind(ints1, ints2), cdij=cdij, verbose)
					
				} else if(identical(ints1, .EmptySet) && !identical(ints2, .EmptySet) && !identical(ints2, .EmptySet) && nr2==1) {
					
					## Case BC2 (1 crossing point in ints2 but none in ints1, bij is within or (aside from the tangent point) entirely outside the window)
					
					wmat[i,j] <- Kisoengine(xi3=xi3, xj3=xj3, win=win2, ints=ints2, cdij=cdij, verbose)

				} else if(identical(ints2, .EmptySet) && !identical(ints1, .EmptySet) && !identical(ints1, .EmptySet) && nr1==1) {
					
					## Case BC2 (1 crossing point in ints2 but none in ints1, bij is within or (aside from the tangent point) entirely outside the window)
					
					wmat[i,j] <- Kisoengine(xi3=xi3, xj3=xj3, win=win1, ints=ints1, cdij=cdij, verbose)

				} else if(nr1 + nr2 <= 2) {
			
					## Cases BC1 and BC3

					wmat[i,j] <- 1
				}

				if(disc) {

					## Also calculate discretized approximation of the same quantity

					ints.disc.ij <- rot.sphere(cbind(rep(Drad[i,j], times=100), lons), northpole=xi, inverse=TRUE)
					circs.disc[i,j] <- sum(in.W(points=ints.disc.ij, win=win))/100
				}
			}
		}
	}

## In the event we want the discretized estimator calculated, the code below turns the output into a data.frame containing the weights matrices
## (each of which could be provided to compileK the argument \code{weights}) for the actual and discretized estimators of Kiso

	if(disc) {
		attr(wmat, "discrete") <- circs.disc
	}
	if(any(is.na(wmat))) {stop("NA in weight matrix")}
	return(wmat)
}


## Function: Kisowedge
## Calculate: Kiso for the wedge (called by Kiso, rather than directly by the user)
## Arguments: X is a two column matrix giving the spherical coordinates of points in the point pattern
##            win is an object of type sphwin
##            r is the vector of distances at which the estimator is calculated
##            nrX, Dmat are respectively the number of rows in X, and the pairwise distance matrix for X
##            rad is the radius of the sphere (it is an element of sphwin but having it separately below is more efficient)
##            disc is a logical, if it is TRUE then the discretized estimation of the weight matrix is performed


Kisowedge <- function(X, win, r, nrX=nrow(X), Dmat=pairdistsph(X), disc=FALSE, verbose=FALSE, rad=win$rad) {

	## Check to ensure the window is a wedge, extract the dihedral angle and save it, its sine and cosine for future use.  Also define output objects.

	stopifnot(win$type=="wedge")
	lon <- win$param[1]
	londiv <- win$param[2]
	slon <- sround(cround(sin(lon)))
	clon <- sround(cround(cos(lon)))
# ajb: 'slon' and 'clon' defined but not used
	wmat <- matrix(nrow=nrX, ncol=nrX)
	ref1 <- win$ref
	ref2 <- c(pi-ref1[1], pi+ref1[2]) %% (2*pi)
	diag(wmat) <- 1
	output <- cbind(0,0)
# ajb: 'output' defined but not used        
	if(disc) {
		circs.disc <- matrix(nrow=nrX, ncol=nrX)
		diag(circs.disc) <- 1
		lons <- seq(0, 2*pi*(1-(100^-1)), length=100)
	}
	Drad <- Dmat/rad
	CD <- sround(cround(cos(Drad)))
	circcens <- rot.sphere(rbind(c(pi/2, londiv+(pi/2)) %% (2*pi), c(pi/2, lon-londiv-(pi/2)) %% (2*pi)), northpole=ref1, inverse=TRUE)
##	print(circcens)
	c12 <- circcens[1,]
	c22 <- circcens[2,]
	c13 <- convert3(c12)
	c23 <- convert3(c22)

	if(lon < pi) {
		win1 <- sphwin(type="band", param=c(0, pi/2), ref=c12)
		win2 <- sphwin(type="band", param=c(0, pi/2), ref=c22)
	} else {
		win1 <- win
		win2 <- win
	}

	for(i in 1:nrX) {

		## define xi and xj as the appropriate objects (spherical coordinates as numeric and matrix, Cartesian coordinates as matrix)

		  xi <- X[i,]
		  ximat <- matrix(xi, nrow=1, ncol=2, byrow=TRUE)
# ajb: 'ximat' defined but not used.                  
		  xi3 <- convert3(xi)
##		  print(xi)
##		  print(xi3)
			for(j in 1:nrX) {
				if(i !=j) {
					## cat("i=", i, ", j=", j, "\n")
					xj <- X[j,]
					xj3 <- convert3(xj)

					## Get dij and its cosine as objects

					dij <- Dmat[i,j]
# ajb: 'dij' is defined but never used                                
					cdij <- CD[i,j]
					
					## Calculate the centres of the caps for which the semicircles from the bases make up the boundary of the wedge

					## Calculate the intersection points
	
					ints1 <- intcircs(n1=c13, n2=xi3, h1=0, h2=cdij)
					ints2 <- intcircs(n1=c23, n2=xi3, h1=0, h2=cdij)
if(verbose){print(ints1)}
if(verbose){print(ints2)}
					nr1 <- nrow(ints1)
					nr2 <- nrow(ints2)

 					## If both arcs have two intersection points, we use Kisowedgeengine to find wij.  Otherwise wij=1.

					if(((identical(.InfiniteSet, ints1) && nr2==2) || (identical(.InfiniteSet, ints2) && nr1==2)) && lon < pi) {

						## Case 1.3
if(verbose){print(1)}
						wmat[i,j] <- 1/2

					} else if(nr1==1 && nr2==2) {

						## Cases W1.4 and W2.3 (ignoring win2 and focussing on win1)

if(verbose){print(2)}
						if(lon >= pi) {

							## If lon < pi, then all points in ints2 are crossing points and are retained; if lon > pi we need to exclude any points in ints2 that aren't 							## crossing points, using Kisowedgeengine
if(verbose){print(2.1)}
							ints2 <- Kisowedgeengine(xi=xi3, xj=xj3, ints=ints2, nr=nr2, c1=c23, c2=c13, ref=ref2, win=win)
						}

						
						if(!identical(ints2, .EmptySet)) {
if(verbose){print(2.11)}
							## If the points in ints2 are crossing points, then we use Kisoengine to find wij.

							wmat[i,j] <- Kisoengine(xi3=xi3, xj3=xj3, ints=ints2, win=win2, cdij=cdij, verbose)
						} else {
if(verbose){print(2.12)}
							## If the points in ints2 are not crossing points, then w_ij=-1

							wmat[i,j] <- 1
						}
					} else if(nr2==1 && nr1==2) {

						## Cases W1.4 and W2.3 (ignoring win1 and focussing on win2)


if(verbose){print(3)}
						if(lon >= pi) {

							## If lon < pi, then all points in ints1 are crossing points and are retained; if lon > pi we need to exclude any points in ints1 that aren't
 							## crossing points, using Kisowedgeengine
if(verbose){print(3.1)}
							ints1 <- Kisowedgeengine(xi=xi3, xj=xj3, ints=ints1, nr=nr1, c1=c13, c2=c23, ref=ref1, win=win, verbose)
						}
						if(!identical(ints1, .EmptySet)) {
							
							## If the points in ints1 are crossing points, then we use Kisoengine to find wij.
if(verbose){print(3.11)}
							wmat[i,j] <- Kisoengine(xi3=xi3, xj3=xj3, ints=ints1, win=win1, cdij=cdij, verbose)
						} else {

							## If the points in ints2 are not crossing points, then w_ij=-1
if(verbose){print(3.12)}
							wmat[i,j] <- 1
						}
					} else if(nr1==2 && nr2==2) {

						## Case W1.5, W2.4
if(verbose){print(4)}
						keptints <- .EmptySet
						keptints1 <- Kisowedgeengine(xi=xi3, xj=xj3, ints=ints1, nr=nr1, c1=c13, c2=c23, ref=ref1, win=win, verbose)
if(verbose){print(keptints1)}
						## If any intersection points in ints2 are crossing point, they are retained in keptints

						if(!identical(keptints1, .EmptySet)) {
if(verbose){print(4.11)}
							keptints <- keptints1
						}
						keptints2 <- Kisowedgeengine(xi=xi3, xj=xj3, ints=ints2, nr=nr2, c1=c23, c2=c13, ref=ref2, win=win, verbose)
if(verbose){print(keptints2)}
						if(!identical(keptints2, .EmptySet)) {

						## If any intersection points in ints2 are crossing point, they are retained in keptints.
if(verbose){print(4.12)}
							if(!identical(keptints, .EmptySet)) {
if(verbose){print(4.121)}
								keptints <- rbind(keptints, keptints2)
							} else {
if(verbose){print(4.122)}
								keptints <- keptints2
							}
if(verbose){print(keptints)}
						}
						if(!identical(keptints, .EmptySet)) {
							
							## If intersection points are retained, use Kisoengine to find wij, otherwise wij=1
if(verbose){print(4.1)}
							wmat[i,j] <- Kisoengine(xi3=xi3, xj3=xj3, ints=keptints, win=win, cdij=cdij, verbose)
						} else {
if(verbose){print(4.2)}
							wmat[i,j] <- 1
						}
					} else if((identical(.EmptySet, ints1) && nr2 == 1) || (identical(.EmptySet, ints2) && nr1 == 1) || ((identical(.InfiniteSet, ints1) || identical(.InfiniteSet, ints2)) && lon >= pi)) {
if(verbose){print(5)}					
					## Cases W1.1, W1.2, W2.1, W2.2
        
						wmat[i,j] <- 1
					} else if(nr1 == 1 && nr2 == 1) {
						wmat[i,j] <- Kisoengine(xi3=xi3, xj3=xj3, ints=rbind(ints1, ints2), win=win, cdij=cdij, verbose)
if(verbose){print(6)}
					}
					## If we want the discretized estimator, we calculate it now.

					if(disc) {
						ints.disc.ij <- rot.sphere(cbind(rep(Drad[i,j], times=100), lons), northpole=xi, inverse=TRUE)
						circs.disc[i,j] <- sum(in.W(points=ints.disc.ij, win=win))/100
					}
				}
			}
	}

	## In the event we want the discretized estimator calculated, the code below turns the output into a data.frame containing the weights matrices
	## (each of which could be provided to compileK the argument ##{weights}) for the actual and discretized estimators of Kiso
	if(disc) {
		attr(wmat, "discrete") <- circs.disc
	}

	wmat
}


## Function: Kisowedgeengine
## Takes crossing points in and determines whether they are crossing points (i.e. in the intersection of the boundary of bij, and the boundary of W (which is a wedge)
## Arguments: xi, xj are 3D vectors (xi centre of bij, xj on circumference of bij)
##	      ints is a 3 column matrix of crossing points
##	      nr is the number of rows in ints
##	      c1, c2 are the centres of the caps that wedge is either the union or intersection of.  c1 is the centre of the cap we are focussed on
##	      ref is a vertex of the wedge
##	      win is the sphwin object defining the wedge
##	      verbose is a logical, if TRUE then output is printed that can help identify bugs

Kisowedgeengine <- function(xi, xj, ints, nr, c1, c2, ref, win, verbose=FALSE) {

	## Create the output object, and orthonormalised vectors
	intsfinal <- matrix(ncol=3)
	e2a <- c2-dot(c1, c2)*c1
	e2 <- e2a/sqrt(sum(e2a^2))
	e3 <- cross(c1, e2)

	## We now find the endpoints of the arc that defines the wedge

	alpha <- atan2(dot(c2,e3), dot(c2,e2)) %% (2*pi)
	phi1 <- (alpha - pi/2) %% (2*pi)
	phi2 <- (alpha +  pi/2) %% (2*pi)
	if(win$param[1] > pi) {
		phi3 <- phi1
		phi1 <- phi2
		phi2 <- phi3
	}
	ref2 <- c(pi-ref[1], pi+ref[2]) %% (2*pi)
if(verbose) {
	print(alpha)
	print(ref2)
	print(phi1)
	print(phi2)
}
	
	## We now determine find the angular coordinate for our intersection points.  The first if loop restricts our focus to intersection points that are not vertices of the wedge and are in W.  ## The second if loop restricts our focus to only those intersection points on the boundary of W.  We retain intersection points that meet these criteria, putting them in a matrix, as Cartesian ## coordinates, and returned.  Kisowedge can give this output to Kisoengine to find wij.

	inWints <- in.W(ints, win)
	for(i in 1:nr) {
		if(!identical(ints[i,], ref) && !identical(ints[i,], ref2) && inWints[i]) {
			phi0 <- atan2(dot(ints[i,], e3), dot(ints[i,], e2)) %% (2*pi)
if(verbose){
 print(phi0)
 print(phi1 < phi2)
 print(phi0 >= phi1)
 print(phi0 <= phi2)
 print(phi1 > phi2)
 print(phi0 <= phi1)
 print(phi0 >= phi2)
}
			if((phi1 < phi2 && phi0 >= phi1 && phi0 <= phi2) || (phi1 >= phi2 && !(phi0 <= phi1 && phi0 >= phi2))) {
				intsfinal <- rbind(intsfinal, ints[i,])
			}
		}
	}
	nri <- nrow(intsfinal)
	if(nri==1) {
		intsfinal <- .EmptySet
	} else {
		intsfinal <- intsfinal[2:nri,]
	}
	return(intsfinal)
}

## Function: Kisopoly
## Calculate: Kiso for the polygon (called by Kiso, rather than directly by the user)
## Arguments: X is a two column matrix giving the spherical coordinates of points in the point pattern
##            win is an object of type sphwin
##            r is the vector of distances at which the estimator is calculated
##            nrX, Dmat are respectively the number of rows in X, and the pairwise distance matrix for X
##            rad is the radius of the sphere (it is an element of sphwin but having it separately below is more efficient)
##            disc is a logical, if it is TRUE then the discretized estimation of the weight matrix is performed



Kisopoly <- function(X, win, r, nrX=nrow(X), Dmat=pairdistsph(X), disc=FALSE, quadwin=NULL, rad=win$rad) {

	## Checking that win is actually a band, and making the bounding circles of colatitude objects

	stopifnot(win$type=="polygon")

	## If we want the discretized estimator, this calls that object	

	if(disc) {
			circs.disc <- matrix(nrow=nrX, ncol=nrX)
			diag(circs.disc) <- 1
			lons <- seq(0, 2*pi*(1-(100^-1)), length=100)
	}

	## Creates the output object

	wmat <- matrix(nrow=nrX, ncol=nrX)
	diag(wmat) <- 1

	## Create other reqired objects

	ref2 <- win$ref2
	cref2 <- sround(cos(ref2))
	sref22 <- (sround(sin(ref2)))^2
# ajb: 'sref22' defined but not used        
	param <- win$param
	nv <- nrow(param)-1
	param1 <- param[c(1:nv,1:2),]
# ajb: 'param1' defined but not used        
	Drad <- Dmat/rad
	CD <- sround(cround(cos(Drad)))
	nv <- nrow(param)-1
	
	## Calculate and store cross products (actual and standardised) of all consecutive pairs of vertices
	
	param3 <- convert3(param)
	param3a <- param3[c(2:nv, 1:2),]
	cvmat <- cbind(param3[,2]*param3a[,3], param3[,3]*param3a[,1], param3[,1]*param3a[,2]) - cbind(param3a[,2]*param3[,3], param3a[,3]*param3[,1], param3a[,1]*param3[,2])
	cvmat <- cvmat[c(1:nv, 1),]
	cvmatn <- cvmat/sqrt(rowSums(cvmat^2))

	for(i in 1:nrX) {

		## Get data for xi in both spherical and Cartesian coordinates

		xi <- X[i,]
		xi3 <- convert3(xi)
		for(j in 1:nrX) {
			if(i != j) {

				## Get data for xj in both spherical and Cartesian coordinates

				xj <- X[j,]
				xj3 <- convert3(xj)

				## Get data for the distance between and dot product of xi and xj


				dij <- Dmat[i,j]
# ajb: 'dij' is defined but never used                                
				cdij <- CD[i,j]
				bijints <- matrix(c(0,0,1), nrow=1, ncol=3, byrow=TRUE)
				for(k in 1:nv) {
					vk03 <- param3[k,]
					vk13 <- param3[k+1,]
					vk23 <- param3[k+2,]
# ajb: 'vk23' defined but not used                                        
					cvn <- cvmatn[k,]	

					## find all crossing points

					ints <- intcircs(n1=cvn, n2=xi3, h1=0, h2=cdij)
					nr <- nrow(ints)
					win1 <- sphwin(type="band", param=c(0, pi/2), ref=convert2(cvn))

					## If there are two crossing points, then we need to use Kisopolyengine to check whether or not they're in W, and if so retain them. (Case P3) 

					if(nr==2) {
						cv2n <- cvmatn[k+1,]
						retpoints <- Kisopolyengine(xi3=xi3, xj3=xj3, vk03=vk03, vk13=vk13, win=win1, ints=ints, cdij=cdij, cvn=cvn, cv2n=cv2n, nr=nr)
						if(!identical(retpoints, .EmptySet)) {
							bijints <- rbind(bijints, retpoints)
						}

					} else if(identical(ints, .InfiniteSet)) {

					## Case P2

					bijints <- rbind(bijints, vk03, vk13)
					} else if(identical(ints, .EmptySet) || nr==1) {

					## Case P1 (no action needed as no points are retained)

					}
				}
				nrb <- nrow(bijints)

				## If the set of retained points (i.e. actual intersection points) is empty, wij=1; otherwise we use Kisoengine to find wij.

				if(nrb==1) {
					wmat[i,j] <- 1
				} else {
					bijints <- bijints[2:nrb,]
					## print(bijints)
					wmat[i,j] <- Kisoengine(xi3=xi3, xj3=xj3, win=win, ints=bijints, cdij=cdij, verbose=FALSE)
					
				}

				if(disc) {

					## Also calculate discretized approximation of the same quantity

					ints.disc.ij <- rot.sphere(cbind(rep(Drad[i,j], times=100), lons), northpole=xi, inverse=TRUE)
					if(is.null(quadwin)) {
						circs.disc[i,j] <- sum(in.W(points=ints.disc.ij, win=win))/100
					} else if(quadwin$type=="quadrangle") {
						circs.disc[i,j] <- sum(in.W(points=ints.disc.ij, win=quadwin))/100
					} else {stop("Internal error: discrete window specification")}
				}
			}
		}
	}

## In the event we want the discretized estimator calculated, the code below turns the output into a data.frame containing the weights matrices
## (each of which could be provided to compileK the argument \code{weights}) for the actual and discretized estimators of Kiso

	if(disc) {
		attr(wmat, "discrete") <- circs.disc
	}
	if(any(is.na(wmat))) {stop("NA in weight matrix")}
	return(wmat)
}


## Function: Kisopolyengine
## Takes crossing points on a circle for which in arc is on the boundary of the polygon W, and determines whether they are crossing points (i.e. in the intersection of the boundary of bij, and the ## boundary of W (which is a polygon)
## Arguments: xi3, xj3 are 3D vectors (xi centre of bij, xj on circumference of bij)
##	      vk03, vk13 are the vertices of the polygon on the circle we are focussed on.  Must be 3D vectors
##	      win is the sphwin object defining the polygon
##	      ints is a 3 column matrix of crossing points
##	      cdij is the dot product of xi3 and xj3
##	      nr is the number of rows in ints
##	      cvn is the centre of the circle we are focussed on (normalized cross product of vk03 and vk13)
##	      cv2n is the centre of the circle which contains the polygon vertex vk13 and the arc between vertex vk13 and the next vertex along from it

Kisopolyengine <- function(xi3, xj3, vk03, vk13, win, ints, cdij=gcdist(xi3, xj3)/win$rad, nr=nrow(ints), cvn=cvn, cv2n=cv2n) {
	intsfinal <- matrix(ncol=3)
	e2a <- cv2n-dot(cvn, cv2n)*cvn
	e2 <- e2a/sqrt(sum(e2a^2))
	e3 <- cross(cvn, e2)
	alpha0 <- atan2(dot(vk03,e3), dot(vk03,e2)) %% (2*pi)
	alpha1 <- atan2(dot(vk13,e3), dot(vk13,e2)) %% (2*pi)
	maxalpha <- max(alpha0, alpha1)
	minalpha <- min(alpha0, alpha1)
	for(i in 1:nr) {
		intsi <- ints[i,]
		alphai <- atan2(dot(intsi,e3), dot(intsi,e2)) %% (2*pi)
		if(maxalpha-minalpha <= pi && maxalpha >= alphai && minalpha <= alphai) {
			intsfinal <- rbind(intsfinal, intsi)
		} else if(((minalpha - maxalpha) %% (2*pi) <=pi) && !(maxalpha > alphai && minalpha < alphai)) {
			intsfinal <- rbind(intsfinal, intsi)	
		}
	}
	nri <- nrow(intsfinal)
	if(nri==1) {
		intsfinal <- .EmptySet
	} else {
		intsfinal <- intsfinal[2:nri,]
	}
	return(intsfinal)
}



## Function: Kisoquad (WORK DISCONTINUED)
## Calculate: Kiso for the quadrangle (called by Kiso, rather than directly by the user)
## Arguments: X is a two column matrix giving the spherical coordinates of points in the point pattern
##            win is an object of type sphwin
##            r is the vector of distances at which the estimator is calculated
##            nrX, Dmat are respectively the number of rows in X, and the pairwise distance matrix for X
##            rad is the radius of the sphere (it is an element of sphwin but having it separately below is more efficient)
##            disc is a logical, if it is TRUE then the discretized estimation of the weight matrix is performed



Kisoquad <- function(X, win, r, nrX=nrow(X), Dmat=pairdistsph(X), disc=FALSE) {

	## Checking that win is a quadrangle, and defining re

	stopifnot(win$type=="quadrangle")
	X <- quadtopoly(X=X, win=win)
	win1 <- X$win
	X <- X$X
	nrX <- nrow(X)
	wmat <- Kisopoly(X=X, win=win1, r=r, nrX=nrX, Dmat=Dmat, disc=disc, quadwin=win)
	return(wmat)
}
	
	## Creating a matrix containing the vertices of the polygon corresponding to the quadrangle

quadtopoly <- function(X=NULL, win=NULL) {
	stopifnot(inherits(X, c("sp2", "sp3")) || ((is.null(X) || inherits(X, "matrix")) && !is.null(win)))
	if(inherits(X, c("sp2", "sp3"))) {
		win <- X$win
		XX <- X$X
	} else if(inherits(X, "matrix") || is.null(X)) {
		XX <- X
	}
	isnX <- ifelse(is.null(X), TRUE, FALSE)
# ajb: 'isnX' is defined but not used.
# ajb: why not just 'isnX <- is.null(X)'        
	rad <- win$rad
	stopifnot(win$type=="quadrangle")
	wp1 <- win$param[1]
	wp2 <- win$param[2]
	wp3 <- win$param[3]
	wp4 <- win$param[4]
	wp43 <- wp4+wp3
	winref1 <- win$ref[1]
	ncX <- ncol(XX)
	vert <- matrix(c(wp1, wp4, wp2, wp4, wp2, wp43, wp1, wp43, wp1, wp4), ncol=2, nrow=5, byrow=TRUE)
	if(wp3 <= pi) {
		winref <- rep(0,4)
	} else {
		winref <- c(0,1,0,1)
	}
	if(sround(winref1*(pi-winref1))!=0) {

		## If the quadvertex (given vertex of the wedge that the quadrangle can be defined as being the union of it and a band) is at neither the north or south pole, we rotate the data to 
		## make the quadvertex the north pole, define the polygon corresponding to this window, then invoke Kisopoly to get the weight matrix wmat

		if(!is.null(XX)) {
			if(ncX==2) {
				XX <- convert3(XX)
				ncX <- 3
			}
			XX <- rot.sphere(XX, c(0,0), inverse=TRUE)
		}
		winref2 <- c(pi/2, wp2, pi/2, wp1)
		winpoly <- sphwin(type="polygon", param=vert, ref=winref, ref2=winref2, rad=rad)
	} else if(sround(winref1-pi)==0) {

		## If the quadvertex is at the south pole, simple changes are required to change the data such that the quadvertex is at the north pole, and to accordingly define the polygon
		## corresponding to this window, then invoke Kisopoly to get the weight matrix wmat

		if(!is.null(XX)) {
			XX[,1] <- pi-XX[,1]
		}
		vert[,2] <- pi-vert[,2]
		winref2 <- c(pi/2, pi-wp2, pi/2, pi-wp1)
		winpoly <- sphwin(type="polygon", param=vert, ref=winref, ref2=winref2, rad=rad)
	} else if(sround(winref1)==0) {

		## If the quadvertex (given vertex of the wedge that the quadrangle can be defined as being the union of it and a band) is at the north pole, we define the polygon corresponding to 
		## this window, then invoke Kisopoly to get the weight matrix wmat
		winref2 <- c(pi/2, wp2, pi/2, wp1)
		winpoly <- sphwin(type="polygon", param=vert, ref=winref, ref2=winref2, rad=rad)
	} else {stop("Internal error: winref1")}

		## If winref1 somehow Dmatoes not satisfy one of the three preceding if arguments, then there must be a bug and the line above ensures that processing Dmatoes not continue
	if(is.null(X)) {
		output <- winpoly
	} else if(inherits(X, "sp2") || (inherits(X, "matrix") && ncX==2)) {
		if(ncX==3) {
			XX <- convert2(XX)
		}
		output <- sp2(X=XX, win=winpoly, check=FALSE)
	} else if(inherits(X, "sp3") || (inherits(X, "matrix") && ncX==3)) {
		output <- sp3(X=XX, win=winpoly, check=FALSE)
	} else stop("Internal error: X")
	return(output)
}


