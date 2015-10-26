nndistsph <-
function(X) {
stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
if(inherits(X, "matrix")) {X <- X} else {
rad <- X$win$rad
X <- X$X}
nrX <- nrow(X)
if(ncol(X)==3) {X <- convert2(X)}
if(nrX>0){
sphdist <- pairdistsph(X)
diag(sphdist) <- Inf
sphmat <- apply(sphdist, 1, min)}
else {
sphmat <- NA}
sphmat
}