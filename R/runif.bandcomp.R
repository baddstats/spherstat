runif.bandcomp <- function(n, win) {
	stopifnot(win$type=="bandcomp")
	win1 <- sphwin(type="band", param=c(0, win$param[1]), ref=win$ref, rad=win$rad)
	win2 <- sphwin(type="band", param=c(win$param[2], pi), ref=win$ref, rad=win$rad)
	aw1 <- area.sphwin(w=win1)
	aw2 <- area.sphwin(w=win2)
	n1 <- rbinom(n=1, size=n, prob=aw1/(aw1+aw2))
	points1 <- runif.band(n=n1, win=win1)
	points2 <- runif.band(n=n-n1, win=win2)
	points <- rbind(points1, points2)
	points
}