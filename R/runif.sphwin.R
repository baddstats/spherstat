runif.sphwin <-
function(n, win=sphwin(type="sphere"), as.sp=TRUE, sp.dim="2") {
stopifnot(inherits(win, "sphwin"))
X <- switch(win$type,
		sphere={runif.sphere(n=n, win=win)},
		band={runif.band(n=n, win=win)},
		bandcomp={runif.bandcomp(n=n, win=win)},
		wedge={runif.wedge(n=n, win=win)},
		polygon={runif.polygon(n=n, win=win)},
		quadrangle={runif.quadrangle(n=n, win=win)},
		)
if(as.sp) {output <- switch(sp.dim,
"2" = {sp2(X=X, win=win, check=FALSE)},
"3" = {sp3(X=convert3(X), win=win, check=FALSE)},
stop("sp.dim must equal 2 or 3")
)} else {output <- X} 
output  
}