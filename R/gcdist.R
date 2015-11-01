gcdist <- function(x,y,rad=1){
  stopifnot(inherits(x, c("matrix", "sp2", "sp3")))
  stopifnot(inherits(y, c("matrix", "sp2", "sp3")))
  if(inherits(x, c("sp2", "sp3"))) {x <- x$X}
  if(inherits(y, c("sp2", "sp3"))) {y <- y$X}
  if(ncol(x) ==2) {x1 <- convert3(x)} else {x1 <- x}
  if(ncol(y) ==2) {y1 <- convert3(y)} else {y1 <- y}
  dp <- x1 %*% t(y1)
  gc <- rad*acos(cround(dp/rad))
  gc
}

gcdistPaired <- function(x, y, rad=1) {
  if(inherits(x, c("sp2", "sp3"))) {x <- x$X} else if(!is.matrix(x))
    stop("x should be a matrix or an object of class sp2 or sp3")
  if(inherits(y, c("sp2", "sp3"))) {y <- y$X} else if(!is.matrix(y))
    stop("y should be a matrix or an object of class sp2 or sp3")
  if(ncol(x) ==2) { x <- convert3(x)} 
  if(ncol(y) ==2) { y <- convert3(y)}
  if(nrow(x) != nrow(y))
    stop("x and y must contain the same numbers of points")
  dp <- rowSums(x * y)
  gc <- rad*acos(cround(dp/rad))
  gc
}

  
