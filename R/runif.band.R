runif.band <-
function(n, win) {
stopifnot(inherits(win, "sphwin") && win$type=="band")
phi <- runif(n=n, min=0, max=2*pi)
theta <- acos(cround(runif(n, cos(win$param[2]), cos(win$param[1]))))
output <- rot.sphere(points=cbind(theta, phi), northpole=win$ref, inverse=TRUE)
output
}
