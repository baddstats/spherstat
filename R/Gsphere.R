Gsphere <- function(X, win=sphwin(type="sphere"),
                    rvals=seq(from=0, to=pi, length=512), ...) {
  stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
  if(inherits(X, "matrix")) {
    stopifnot(inherits(win, "sphwin"))
  } else {
    win <- X$win
    X <- X$X
    stopifnot(inherits(win, "sphwin"))
  }
  rad <- win$rad
  if(length(X)==0) {
    stop("Cannot estimate G function of an empty set")
  } else {
    D <- nndistsph(X=X)
    B <- bdist.sphwin(X=X, win=win)
    lambda <- intensitysph(X=X, win=win)
    r <- rad*rvals
    h <- eroded.areas.sphwin(win=win, r=r, ...)
    if(any(!is.finite(h))) {
      stop("Error: h not finite")
    }
    f <- compileCDF(D=D,B=B, r=rvals, han.denom=h, check=FALSE)
    dn <- fvnames(f, ".")
    f <- rebadge.fv(f, new.fname="Gsphere")
    lr <- length(r)
    n <- length(D)
    nnmat <- matrix(rep(D, lr), ncol=n, nrow=lr, byrow=TRUE)
    frsmod <- (area.sphwin(win)*rowSums(nnmat<=r))/(h*n)
    f <- bind.fv(f,
                 data.frame(rsmodif=frsmod),
                 labl="%s[bordmodif](r)",
                 desc="modified border corrected estimate of %s")
    ftheo <- 1 - exp(-lambda *2 * pi * (rad^2) * (1-cos(r/rad)))
    f <- bind.fv(f, data.frame(theo=ftheo),
                 labl="%s[theo](r)", desc="theoretical value of %s")
    fvnames(f, ".") <- c(dn, "rsmodif", "theo")
  }
  return(f)
}
