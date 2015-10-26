gcdist <-
function(x,y,rad=1){
stopifnot(inherits(x, c("matrix", "sp2", "sp3")) && inherits(y, c("matrix", "sp2", "sp3")))
if(inherits(x, c("sp2", "sp3"))) {x <- x$X}
if(inherits(y, c("sp2", "sp3"))) {y <- y$X}
if(ncol(x) ==2) {x1 <- convert3(x)} else {x1 <- x}
if(ncol(y) ==2) {y1 <- convert3(y)} else {y1 <- y}
dp <- x1 %*% t(y1)
gc <- rad*acos(cround(dp/rad))
gc
}
