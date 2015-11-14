polysph.area.W <- function(win) {
	stopifnot(inherits(win, c("sp2", "sp3", "sphwin")))
	if(inherits(win, "sphwin")) {
		win <- win
	}
	else {
		win <- win$win
	}
	stopifnot(win$type=="polygon")
	rad <- win$rad
	angs <- sph.angles(win=win)
	a.polygon <- (rad^2)*(sum(angs) - (length(angs)-2)*pi)
	a.polygon
}