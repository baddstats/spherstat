Gsphere <- function(X, win=sphwin(type="sphere"),
			r=NULL, ...) {
	stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
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
	if(nrow(X)!=0) {
		D <- nndistsph(X=X)
		B <- bdist.sphwin(X=X, win=win)
		lambda <- intensitysph(X=X, win=win)
		h <- eroded.areas.sphwin(win=win, r=r, ...)
		if(any(!is.finite(h))) {
			stop("Error: eroded areas function not finite")
		}
		f <- compileCDF(D=D,B=B, r=r, han.denom=h, check=FALSE)
		dn <- fvnames(f, ".")
		f <- rebadge.fv(f, new.fname="G")
		lr <- length(r)
		n <- length(D)
		nnmat <- matrix(rep(D, lr), ncol=n, nrow=lr, byrow=TRUE)
		frsmod <- (area.sphwin(win)*rowSums(nnmat<=r))/(h*n)
		f <- bind.fv(f,
			data.frame(rs.modif=frsmod),
			labl="%s[bord.modif](r)",
			desc="modified border corrected estimate of %s")
		ftheohaz <- 2*pi*lambda*(rad^2)*sin(r/rad)
		f <- bind.fv(f, data.frame(theohaz=ftheohaz), labl="lambda[pois](r)", "theoretical Poisson hazard h(r)")
		fvnames(f, ".") <- c(dn, "theohaz")
		ftheo <- 1 - exp(-lambda *2 * pi * (rad^2) * (1-cos(r/rad)))
		f <- bind.fv(f, data.frame(theo=ftheo),
			labl="%s[theo](r)", desc="theoretical value of %s")
		fvnames(f, ".") <- c(dn, "rs.modif", "theo")
	}
	else {
		Gs <- rep(0, 512)
		df <- data.frame(r=r, km=Gs, hazard=Gs, han=Gs, rs=Gs, rs.modif=Gs, raw=Gs, theohaz=Gs, theo=Gs)
		alim <- range(df$r[df[["km"]] <= 0.9])
		nama <- c("r", "km", "hazard", "han", "rs", "rs.modif", "raw", "theohaz", "theo")
		iscdf <- c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)
		labl <- c("r", "hat(%s)[km](r)", "lambda(r)", "hat(%s)[han](r)", 
		"hat(%s)[bord](r)", "hat(%s)[bord.modif](r)", "hat(%s)[raw](r)", "lambda[pois](r)", "hat(%s)[theo](r)")
		desc <- c("distance argument r", "Kaplan-Meier estimate of %s", 
		"Kaplan-Meier estimate of hazard function lambda(r)", 
		"Hanisch estimate of %s", "border corrected estimate of %s", "modified border corrected estimate of %s",
		"uncorrected estimate of %s", "theoretical Poisson hazard h(r)", "theoretical value of %s")
		f <- fv(df, "r", substitute(G(r), NULL), "km", . ~ r, 
			alim, labl, desc, fname = "G")
		fvnames(f, ".") <- nama[iscdf]
	}
	return(f)
}
