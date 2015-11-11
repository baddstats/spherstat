sphsin <- function(d1, d2, theta1, theta2, rad=1) {
	stopifnot((is.null(d2) && !is.null(theta2)) || (!is.null(d2) && is.null(theta2)))
	if(is.null(d2)) {
		stopifnot(sround(sin(theta1)) != 0)
		ss <- cround(sin(d1/rad)*sin(theta2)/sin(theta1))
		output <- rad*asin(ss)
	} else if(is.null(theta2)) {
		stopifnot(sround(sin(d1)) != 0)
		ss <- cround(sin(d2/rad)*sin(theta1)/sin(d1))
		output <- asin(ss)
	}
	output
}
