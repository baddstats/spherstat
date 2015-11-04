Fsphere <- function(X, refpoints,
                    win=sphwin(type="sphere"),
                    r=NULL, ...) {
  stopifnot(inherits(X, c("sp2", "sp3", "matrix")) &&
            inherits(refpoints, c("sp2", "sp3", "matrix")))
  if(inherits(X, "matrix")) {
    stopifnot(inherits(win, "sphwin"))
  } else {
    win <- X$win
    X <- X$X
    stopifnot(inherits(win, "sphwin"))
  }
  rad <- win$rad
  if(is.null(r)) {
    rmax <- rmax.rule.sphwin(win)
    r <- seq(0, rmax, length=512)
  }
  if(inherits(refpoints, c("sp2", "sp3")) &&
     refpoints$win$rad != rad) {
    stop("X and refpoints have different radii")
  }
  if(inherits(refpoints, "matrix")) {
    refpoints <- refpoints
  } else {
    refpoints <- refpoints$X
  }
  D <- nncrosssph(X=X, Y=refpoints)
  B <- bdist.sphwin(X=refpoints, win=win)
  lambda <- intensitysph(X=X, win=win)
  h <- eroded.areas.sphwin(win=win, r=r, ...)
  if(any(!is.finite(h))) {
    stop("Error: h not finite")
  }
  f <- compileCDF(D=D,B=B, r=r, han.denom=h, check=FALSE)
  dn <- fvnames(f, ".")
  f <- rebadge.fv(f, new.fname="Fsphere")
  ftheo <- 1 - exp(-lambda *2 * pi * (rad^2) * (1-cos(r/rad)))
  f <- bind.fv(f, data.frame(theo=ftheo),
               labl="%s[theo](r)",
               desc="theoretical value of %s")
  fvnames(f, ".") <- c(dn, "theo")
  return(f)
}
