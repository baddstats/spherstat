convert.globe <- function(X) {
  stopifnot(inherits(X, c("sp2", "sp3", "matrix")))
  if(!inherits(X, "matrix")) 
    X <- X$X
  if(ncol(X)==3) {
    X <- convert2(X)
  }
  Xdegrees <- X*(180/pi)
  Xlat <- 90 - Xdegrees[,1]
  Xlon <- Xdegrees[,2]
  Xlon <- ifelse(Xlon <= 180, Xlon, Xlon - 360)
  Xglobe <- cbind(lon=Xlon, lat=Xlat)
  return(Xglobe)
}
