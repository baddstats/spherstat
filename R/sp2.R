sp2 <- function(X, win=sphwin(type="sphere"), check=TRUE) {
	stopifnot(ncol(X)==2 && inherits(win, "sphwin"))
	if(check) {
		n <- nrow(X)
		stopifnot(sum(in.W(points=X, win=win))==n && sum(X >= 0)==2*n && sum(X[,1] <= pi)==n && sum(X[,2] < 2*pi)==n)
	}
	result <- list(X=X, win=win)
	class(result) <- c("sp2", class(result))
	result
}