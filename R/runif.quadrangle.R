runif.quadrangle <-
function(n, win) {
stopifnot(inherits(win, "sphwin") && win$type == "quadrangle")
phi1 <- runif(n, 0, win$param[3])
phi.rot1 <- (phi1 + win$param[4]) %% (2*pi)
theta <- acos(cround(runif(n, cos(win$param[2]), cos(win$param[1]))))
output <- cbind(theta, phi.rot1)
output <- rot.sphere(points=output, northpole=win$ref, inverse=TRUE)
output
}
