summary.sphwin <- function(object, ...) {
	stopifnot(inherits(object, "sphwin"))
	print.sphwin(object)
	print(cat("\n", "Window area: ", area.sphwin(object), " square units.", "\n"))
}
