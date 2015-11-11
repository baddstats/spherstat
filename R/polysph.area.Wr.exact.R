polysph.area.Wr.exact <- function(win, r=0) {
	stopifnot(inherits(win, c("sp2", "sp3", "sphwin")))
	if(inherits(win, "sphwin")) {
		win <- win
	}
	else {
		win <- win$win
	}
	stopifnot(win$type=="polygon" || win$type == "wedge")
	rad <- win$rad
	lr <- length(r)-1
	areaWr <- switch(win$type,
		wedge = {
			phi <- pi/2 - r[2:lr]/rad
			theta <- win$param[1]/2
			h <- pi - 2*sphsin(d1=r[2:lr], d2=NULL, theta1=theta, theta2=pi/2, rad=rad)
			psi1 <- sphcos(d1=phi, d2=h, d3=phi, theta=NULL, rad=1)
			psi0 <- sphcos(d1=phi, d2=phi, d3=h, theta=NULL, rad=1)
			area.wedge <- (rad^2)*psi0*(1-cos(phi))
			area.tri <- (rad^2)*(psi0 + 2*psi1 - pi)
			areaWr1 <- c(area.sphwin(w=win), 2*(area.wedge-area.tri), 0)
			areaWr1
		},
		stop("Unsupported window type")
	)
	areaWr
}
