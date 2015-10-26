summary.sp2 <- function(object, ...) {
stopifnot(inherits(object, "sp2"))
nrX <- nrow(object$X)
cat("Spherical point pattern: ", nrX, " points.", "\n")
cat("Average intensity ", intensitysph(object), " points per square unit.", "\n")
cat("\n", "Points are given in spherical coordinates i.e. (colatitude, longitude).", "\n")
cat("\n", summary.sphwin(object$win))
}