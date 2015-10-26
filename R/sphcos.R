sphcos <-function(d1, d2, d3, theta, rad=1) {
	stopifnot((is.null(d3) && !is.null(theta)) || (!is.null(d3) && is.null(theta)))
	if(is.null(theta)) {
		sc <- cround((cos(d3/rad)-cos(d1/rad)*cos(d2/rad))/(sin(d1/rad)*sin(d2/rad)))
		output <- acos(sc)
	} else if (is.null(d3)) {
		sc <- cround(cos(d1/rad)*cos(d2/rad) + sin(d1/rad)*sin(d2*rad)*cos(theta))
		output <- rad*acos(sc)
}
output
}
