Fsphere <- function(X, refpoints=NULL,
	win=sphwin(type="sphere"),
	r=NULL, ...) {
	stopifnot(inherits(X, c("sp2", "sp3", "matrix")) &&
			inherits(refpoints, c("sp2", "sp3", "matrix", "NULL")))
	if(inherits(X, "matrix")) {
		stopifnot(inherits(win, "sphwin"))
	}
	else {
		win <- X$win
		X <- X$X
		stopifnot(inherits(win, "sphwin"))
	}
	rad <- win$rad
	if(is.null(r)) {
    	rmax <- rmax.rule.sphwin(win)
	r <- seq(0, rmax, length=512)
	}
	if(is.null(refpoints)) {
		refpoints <- rpoispp.sphwin(lambda=100/area.sphwin(win), win=win)
	}
	if(inherits(refpoints, c("sp2", "sp3")) &&
		refpoints$win$rad != rad) {
 		stop("X and refpoints have different radii")
	}
	if(inherits(refpoints, "matrix")) {
		refpoints <- refpoints
	}
	else {
		refpoints <- refpoints$X
	}
	if(nrow(X) != 0) {
		D <- nncrosssph(X=X, Y=refpoints)
		B <- bdist.sphwin(X=refpoints, win=win)
		lambda <- intensitysph(X=X, win=win)
		h <- eroded.areas.sphwin(win=win, r=r, ...)
		if(any(!is.finite(h))) {
		stop("Error: h not finite")
	}
	f <- compileCDF(D=D,B=B, r=r, han.denom=h, check=FALSE)
	dn <- fvnames(f, ".")
	f <- rebadge.fv(f, new.fname="F")
	ftheohaz <- 2*pi*lambda*(rad^2)*sin(r/rad)
	f <- bind.fv(f, data.frame(theohaz=ftheohaz), labl="lambda[pois](r)", "theoretical Poisson hazard h(r)")
	fvnames(f, ".") <- c(dn, "theohaz")
	ftheo <- 1 - exp(-lambda *2 * pi * (rad^2) * (1-cos(r/rad)))
	f <- bind.fv(f, data.frame(theo=ftheo),
		labl="%s[theo](r)",
		desc="theoretical value of %s")
	fvnames(f, ".") <- c(dn, "theo")
	}
	else {
		Fs <- rep(0, 512)
		df <- data.frame(r=r, km=Fs, hazard=Fs, han=Fs, rs=Fs, raw=Fs, theohaz=Fs, theo=Fs)
		alim <- range(df$r[df[["km"]] <= 0.9])
		nama <- c("r", "km", "hazard", "han", "rs", "raw", "theohaz", "theo")
		iscdf <- c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE)
		labl <- c("r", "hat(%s)[km](r)", "lambda(r)", "hat(%s)[han](r)", 
		"hat(%s)[bord](r)", "hat(%s)[raw](r)", "lambda[pois](r)", "hat(%s)[theo](r)")
		desc <- c("distance argument r", "Kaplan-Meier estimate of %s", 
		"Kaplan-Meier estimate of hazard function lambda(r)", 
		"Hanisch estimate of %s", "border corrected estimate of %s", 
		"uncorrected estimate of %s", "theoretical Poisson hazard h(r)", "theoretical value of %s")
		f <- fv(df, "r", substitute(F(r), NULL), "km", . ~ r, 
			alim, labl, desc, fname = "F")
		fvnames(f, ".") <- nama[iscdf]

	}
	return(f)
}
