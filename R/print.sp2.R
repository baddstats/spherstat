print.sp2 <- function(x, ...) {
	stopifnot(inherits(x, "sp2"))
	nrX <- nrow(x$X)
	cat("Spherical point pattern: ", nrX, " points.", fill=TRUE)
	cat("\n", print.sphwin(x$win))
}
