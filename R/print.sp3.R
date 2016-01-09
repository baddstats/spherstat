print.sp3 <- function(x, ...) {
	stopifnot(inherits(x, "sp3"))
	nrX <- nrow(x$X)
	cat("Spherical point pattern: ", nrX, " points.", fill=TRUE)
	cat("\n", print.sphwin(x$win))
}
