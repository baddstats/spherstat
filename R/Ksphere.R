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


Ksphere <- function(X, win=sphwin(), rvals,
                    correction=c("un", "iso", "rs", "rsm"),
                    ratio=FALSE) {

  ## standard checks of the X and win arguments
  stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
  if(!inherits(X, "matrix")) {
    win <- X$win
    X <- X$X
  }
  stopifnot(inherits(win, "sphwin"))

  ## Create objects that will be used in later calculaations
  
  r <- rvals
  rad <- win$rad
  nrX <- nrow(X)
  lr <- length(r)
  areaX <- area.sphwin(win)
  lambda <- nrX/areaX
  lambda2 <- nrX*(nrX-1)/areaX^2
  denom <- lambda2 * areaX
  Dmat <- pairdistsph(X)

  ## Create fv object that contains the theoretical value of K

  Kdf <- data.frame(r = r, theo = 2*pi *rad*(1-cos(r/rad)))
  desc <- c("distance argument r", "theoretical Poisson %s")
  K <- ratfv(Kdf, NULL, denom, "r", quote(K(r)), "theo",
             . ~ r, range(r), c("r", "%s[pois](r)"),
             desc, fname = "K", ratio = ratio)

  ## The un, rs and rsm estimates require similar calculations
  ## so we have an if loop covering them all, 

  if(any(correction=="un") || any(correction=="rs") || any(correction=="rsm")) {

    ## We create a matrix xdij, for which the [j,i]th cell of xdij
    ## gives the number of points within distance r[i] of point X[j]

    xdij <- matrix(ncol=lr, nrow=nrX)
    for(i in 1:lr) {
      xdij[,i] <- colSums(Dmat <= r[i]) - 1
    }
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

      ## We create rmat, which in each row lists all value of r.
      ## Then we create bX, a matrix where all the values in the jth row are
      ## the minimum distance from X[j] to the boundary of the window.

      rmat <- matrix(rep(r, each=nrX), nrow=nrX, ncol=lr, byrow=FALSE)
      bX <- matrix(rep(bdist.sphwin(X=X, win=win), each=lr),
                   nrow=nrX, ncol=lr, byrow=TRUE)

      ## xinWr is a matrix where the [j,i]th cell takes value 1
      ## if X[j] is in W_-r[i] and 0 otherwise.
      ## The ith value nXwr is the number of points
      ## in W_-r[i].
      ## Kbnum a matrix where the [j,i]th cell takes value 0
      ## if X[j] is not in W_-r[i], and the number of points with distance r[i]
      ## of X[j] otherwise.
      ## Kbcolsum is the list of values of the numerator of the border
      ## and modified border corrected estimates of K for all r.

      xinWr <- (bX >= rmat)
      nXwr <- colSums(xinWr)
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

    Kiso <- Kiso(X=X, win=win, r=r, Dmat=Dmat, nrX=nrX, denom=denom)
    Kiso <- data.frame(r=r, iso=Kiso$est)
    Kiso <- fv(Kiso, "r", quote(hat(K(r))), "iso", NULL,
               range(r), c("r", "%s[iso](r)"),
               desc=c("distance argument r",
                 "isotropic-corrected estimate of %s"), fname="K")
    K <- bind.fv(K, Kiso)
    fvnames(K, ".y") <- "iso"
  }
  return(K)
}
