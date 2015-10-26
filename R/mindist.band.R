mindist.band <-
function(X, win) {
stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
if(inherits(X, "matrix")) {stopifnot(inherits(win, "sphwin"))} else {
win <- X$win
X <- X$X
}
stopifnot(win$type=="band" || win$type=="bandcomp")
rad <- win$rad
n <- nrow(X)
X1 <- matrix(nrow=n, ncol=2)
X1[,1] <- abs(X[,1]-win$param[1])
X1[,2] <- abs(win$param[2]-X[,1])
mindists <- apply(X1, 1, min)*rad
mindists
}
