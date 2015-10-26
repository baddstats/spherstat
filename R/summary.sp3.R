summary.sp3 <- function(object, ...) {
stopifnot(inherits(object, "sp3"))
nrX <- nrow(object$X)
cat("Spherical point pattern: ", nrX, " points.", "\n")
cat("Average intensity ", intensitysph(object), " points per square unit.", "\n")
cat("\n", "Points are given in Cartesian coordinates.", "\n")
cat("\n", summary.sphwin(object$win))
}