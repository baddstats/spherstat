sp3 <- function(X=X, win=sphwin(type="sphere"), check=TRUE) {
stopifnot(ncol(X)==3 && inherits(win,"sphwin"))
if(check) {
stopifnot(in.W(points=X, win=win))
n <- nrow(X)
pointrad <- (X[,1]^2)+(X[,2]^2)+(X[,3]^2)
rad1 <- rep(win$rad, n)
stopifnot(sum(abs(rad1-pointrad)) <= (10^-16)*n)
}
result <- list(X=X, win=win)
class(result) <- c("sp3", class(result))
result
}