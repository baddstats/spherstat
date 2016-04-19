
## File Ksphere.R
## (was: '20140723 Ksphere.txt')

## Function: Ksphere
## Calculates: estimates of the K function for a point pattern
## Arguments: X is the point pattern.  Can be a matrix, or object of class sp2 or sp3
##	      win is the window.  Not required and ignored if X is an object of class sp2 or sp3
##	      r is a numeric containing the distances at which the estimate(s) of K are to be calculated
##	      correction lists the types of estimates to be calculated.  items can include "un" (uncorrected/raw estimate of K), 
##	      "iso" (isotropic-corrected estimate of K), "rs" (border corrected or reduced sample estimator of K) and 
##	      "rsm" (modified border corrected or reduced sample estimator of K)
##	      ratio is an internal argument.
## Output: an fv object containing the chosen estimates of K for X.


Ksphere <- function(X, win=sphwin(), r=NULL,
                    correction=c("un", "iso", "rs", "rsm"),
                    ratio=FALSE, lambdavalues=NULL, update=TRUE) {

	## standard checks of the X and win arguments
	stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
	if(!inherits(X, "matrix")) {
		win <- X$win
		X <- X$X
	}
	stopifnot(inherits(win, "sphwin"))

	if(is.null(r)) {
	rmax <- rmax.rule.sphwin(win)
	r <- seq(0, rmax, length=512)
	}

	correction <- match.arg(correction,
			c("un", "iso", "rs", "rsm", "best"),
			several.ok=TRUE)
	correction[correction == "best"] <- "iso"

	if(!is.null(lambdavalues) && !all(correction == "iso")) {
		warning("For inhomogeneous K, only the isotropic correction is available")
		correction <- "iso"
	}

	if(!is.null(lambdavalues) && inherits(lambdavalues, "sphppm")) {
		model <- lambdavalues
		if(update) {
			## refit model to current data
			XX <- sp2(X, win)
			model <- update(model, XX)
		}
		lambdavalues <- fitted(model)
	}

	## Create objects that will be used in later calculations
	rad <- win$rad
	nrX <- nrow(X)
	lr <- length(r)
	areaX <- area.sphwin(win)
	if(is.null(lambdavalues)) {
		lambda <- nrX/areaX
		lambda2 <- nrX*(nrX-1)/areaX^2
		denom <- lambda2 * areaX
	}
	else {
		denom <- areaX
	}

	Dmat <-	if(all(correction == "iso") && win$type %in% c("band", "bandcomp")) 
			NULL else pairdistsph(X)
	## Create fv object that contains the theoretical value of K
	Kdf <- data.frame(r = r, theo = 2*pi *rad*(1-cos(r/rad)))
	desc <- c("distance argument r", "theoretical Poisson %s")
	K <- ratfv(Kdf, NULL, denom, "r", quote(K(r)), "theo",
	 ~ r, range(r), c("r", "%s[pois](r)"),
	desc, fname = "K", ratio = ratio)
	
	## We now give the code for the cases where X contains points, and where X is empty

	if(nrow(X) !=0) {
		## The un, rs and rsm estimates require similar calculations
		## so we have an if loop covering them all, 

		if(any(correction=="un") || any(correction=="rs") || any(correction=="rsm")) {

			## We create a matrix xdij, for which the [j,i]th cell of xdij
			## gives the number of points within distance r[i] of point X[j]
			xdij <- neighbourcount(Dmat, r)	

			## WAS:
			## xdij <- matrix(ncol=lr, nrow=nrX)
			## for(i in 1:lr) {
			##  xdij[,i] <- colSums(Dmat <= r[i]) - 1
			## }

	
			if(any(correction=="un")) {
				
				## The uncorrected estimator is simply xdij divided by the
				## squared intensity of the point pattern
				Kraw <- data.frame(r = r, un = colSums(xdij)*areaX/(nrX*(nrX-1)))
				Kraw <- fv(Kraw, "r", quote(hat(K(r))), "un",
					NULL, range(r), c("r", "%s[un](r)"),
					desc=c("distance argument r", "uncorrected estimate of %s"),
					fname="K")
				K <- bind.fv(K, Kraw)
				fvnames(K, ".y") <- "un"
			}	

			if(any(correction=="rs") || any(correction=="rsm")) {	

				## xinWr is a matrix where the [i,j]th cell takes value TRUE
				## if X[i] is in W_(-r[j]) and FALSE otherwise.
				bdX <- bdist.sphwin(X=X, win=win)
				xinWr <- outer(bdX, r, ">=")
      
				## The ith value nXwr is the number of points
				## in W_(-r[i]).
				nXwr <- colSums(xinWr)

				## Kbnum a matrix where the [j,i]th cell takes value 0
				## if X[j] is not in W_-r[i], and the number of points with distance r[i]
				## of X[j] otherwise.
				## Kbcolsum is the list of values of the numerator of the border
				## and modified border corrected estimates of K for all r.
				Kbnum <- xinWr * xdij
				Kbcolsum <- colSums(Kbnum)
				if(any(correction=="rs")) {
					
					## We find the denominator of the border-corrected estimator,
					## and calculate K, using 0/0 = 0 to keep all values of the estimator
					## finite.  We then add this to the output object	

					Krs <- Kbcolsum*areaX/(nrX*nXwr)
					Krs[nXwr==0] <- 0
					Krs <- data.frame(r=r, border=Krs)
					Krs <- fv(Krs, "r", quote(hat(K(r))), "border", NULL,
						range(r), c("r", "%s[bord](r)"),
						desc=c("distance argument r",
						"border-corrected estimate of %s"), fname="K")
					K <- bind.fv(K, Krs)
					fvnames(K, ".y") <- "border"
				}
				if(any(correction=="rsm")) {	

					## Same drill as for the border-corrected estimator,
					## but for the modified border corrected estimator	

					Krsm <- Kbcolsum*(areaX^2)/((nrX^2)*eroded.areas.sphwin(win, r))
					Krsm[1] <- 0
					Krsm <- data.frame(r=r, bord.modif=Krsm)
					Krsm <- fv(Krsm, "r", quote(hat(K(r))), "bord.modif", NULL,
						range(r), c("r", "%s[bord.modif](r)"),
						desc=c("distance argument r",
						"modified border-corrected estimate of %s"), fname="K")
					K <- bind.fv(K, Krsm)
					fvnames(K, ".y") <- "bord.modif"
				}
			}
		}
		if(any(correction=="iso")) {

			## If the isotropic corrected estimator is wanted,
			## this code invokes the relevant code.

			Kiso <- Kiso(X=X, win=win, r=r, Dmat=Dmat, nrX=nrX, denom=denom,
				lambda=lambdavalues)
			Kiso <- data.frame(r=r, iso=Kiso$est)
			Kiso <- fv(Kiso, "r", quote(hat(K(r))), "iso", NULL,
				range(r), c("r", "%s[iso](r)"),
				desc=c("distance argument r",
				"isotropic-corrected estimate of %s"), fname="K")
			K <- bind.fv(K, Kiso)
			fvnames(K, ".y") <- "iso"
		}
	}
	else {
		if(any(correction=="un")) {
			Kraw <- data.frame(r = r, un = rep(NaN, 512))
			Kraw <- fv(Kraw, "r", quote(hat(K(r))), "un",
				NULL, range(r), c("r", "%s[un](r)"),
				desc=c("distance argument r", "uncorrected estimate of %s"),
				fname="K")
			K <- bind.fv(K, Kraw)
			fvnames(K, ".y") <- "un"
		}
		if(any(correction=="rs")) {
			Krs <- data.frame(r = r, border = rep(NaN, 512))
			Krs <- fv(Krs, "r", quote(hat(K(r))), "border", NULL,
				range(r), c("r", "%s[bord](r)"),
				desc=c("distance argument r",
				"border-corrected estimate of %s"), fname="K")
			K <- bind.fv(K, Krs)
			fvnames(K, ".y") <- "border"
		}
		if(any(correction=="rsm")) {
			Krsm <- data.frame(r = r, bord.modif = rep(NaN, 512))
			Krsm <- fv(Krsm, "r", quote(hat(K(r))), "bord.modif", NULL,
				range(r), c("r", "%s[bord.modif](r)"),
				desc=c("distance argument r",
				"modified border-corrected estimate of %s"), fname="K")
			K <- bind.fv(K, Krsm)
			fvnames(K, ".y") <- "bord.modif"
		}
		if(any(correction=="iso")) {
			Kiso <- data.frame(r = r, iso = rep(NaN, 512))
			Kiso <- fv(Kiso, "r", quote(hat(K(r))), "iso", NULL,
				range(r), c("r", "%s[iso](r)"),
				desc=c("distance argument r",
				"isotropic-corrected estimate of %s"), fname="K")
			K <- bind.fv(K, Kiso)
			fvnames(K, ".y") <- "iso"
		}
	}

	return(K)
}


neighbourcount <- function(Dmat, r) {
	## Dmat is a matrix of pairwise distances
	nD <- nrow(Dmat)
	## don't count identical pairs
	diag(Dmat) <- Inf
	## r is a vector of distance thresholds
	nr <- length(r)
	## For each distance value, find the r interval (r[k-1], r[k]] containing it
	rindex <- nr-findInterval(-Dmat, -rev(r))
	rindex <- matrix(rindex, nrow=nD, ncol=nD)
	## Cumulative count of distance values in each r interval 
	nc <- apply(rindex, 1, cumtab, m=0:(nr-1))
	## Return a matrix whose [i,j] entry equals the number
	## of points lying within distance r[j] of point i.
	return(t(nc))
}

cumtab <- function(z, m) { cumsum(table(factor(z, levels=m))) }
Lsphere <- function(X, ...) {
	K <- Ksphere(X, ...)
	L <- eval.fv(sqrt(K))
	L <- rebadge.fv(L, quote(L(r)), "L", names(K), new.labl = attr(K, 
        "labl"))
	return(L)
}


Kpp <- function(X, ...) {
	K <- Ksphere(X, ...)
	Kpp <- eval.fv(acos(1-K/(2*pi)))
	Kpp <- rebadge.fv(Kpp, quote(K_pp(r)), "K_pp", names(K), new.labl = attr(K, 
        "labl"))
	return(Kpp)
}
