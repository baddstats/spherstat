rFisher <- function(n, mode, kappa, win=sphwin()) {
stopifnot(inherits(win, "sphwin") && kappa > 0)
lambda <- exp(-2*kappa)
r1 <- runif(n,0,1)
r2 <- runif(n,0,1)
theta <- 2*asin(sqrt(-log(lambda + (1-lambda)*r1)/(2*kappa)))
theta <- theta
phi <- 2*pi*r2
rf1 <- cbind(theta, phi)
rf1 <- rot.sphere(points=rf1, northpole=mode, inverse=TRUE)
rf2 <- sp2(X=rf1, win=win)
rf2
}
