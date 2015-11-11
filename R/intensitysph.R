intensitysph <- function(X, win=sphwin(type="sphere"), ...) {
	stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
	if(inherits(X, "matrix")) {
		stopifnot(inherits(win, "sphwin"))
	}
	else {
		win <- X$win
		X <- X$X
		stopifnot(inherits(win, "sphwin"))
	}
	lambda <- nrow(X)/area.sphwin(win)
	lambda
}
