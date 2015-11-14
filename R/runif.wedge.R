runif.wedge <- function(n, win) {
	stopifnot(inherits(win, "sphwin") && win$type == "wedge")
	theta <- acos(runif(n, -1, 1))
	phi1 <- runif(n, 0, win$param[1])
	phi.rot1 <- (phi1 + win$param[2]) %% (2*pi)
	output <- rot.sphere(cbind(theta, phi.rot1), northpole=win$ref, inverse=TRUE)
	output
}
