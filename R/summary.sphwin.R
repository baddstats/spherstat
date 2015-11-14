summary.sphwin <- function(object, ...) {
	stopifnot(inherits(object, "sphwin"))
	print.sphwin(object)
	cat("\n", "Window area: ", area.sphwin(object), " square units.", "\n")
}
