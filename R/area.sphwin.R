area.sphwin <- function(w) {
	stopifnot(inherits(w, c("sp2", "sp3", "sphwin")))
	if(inherits(w, "sphwin")) {
		w <- w
	}
	else {
		w <- w$w
	}
	rad <- w$rad
	a.sph <- switch(w$type,
			sphere = {
				4*pi*(rad^2)
			},
			band = {
				2*(rad^2)*pi*(cos(w$param[1])-cos(w$param[2]))
			},
			bandcomp = {
				2*(rad^2)*pi*(2-cos(w$param[1])+cos(w$param[2]))
			},
			wedge = {
				w$param[1]*2*(rad^2)
			},
			polygon = {
				polysph.area.W(win=w)
			},
			quadrangle = {
				(rad^2)*w$param[3]*(cos(w$param[1])-cos(w$param[2]))
			},
			{
			stop("Unrecognised window type")
			}
		)
	a.sph
}
