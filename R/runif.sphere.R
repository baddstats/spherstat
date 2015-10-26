runif.sphere <-
function(n, win=sphwin(type="sphere")) {
stopifnot(inherits(win, "sphwin") && win$type == "sphere")
theta <- acos(runif(n, -1, 1))
phi <- runif(n, 0, 2*pi)
output <- cbind(theta, phi)
output
}
