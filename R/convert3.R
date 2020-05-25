convert3 <- function(points, rad=1) {
  stopifnot(length(points) == 2 || ncol(points)==2 || inherits(points, "sp2"))
  if(inherits(points, "sp2")) {
    rad <- points$win$rad
    points <- points$X
  }
  if(!is.matrix(points)) {
    points <- matrix(points, ncol=2, byrow=TRUE)
  } 
  stopifnot(ncol(points)==2)
  theta <- points[,1]
  phi <- points[,2]
  sintheta <- sin(theta)
  points3 <-  cbind(sintheta * cos(phi),
                    sintheta * sin(phi),
                    cos(theta))
  if(rad != 1) 
        points3 <- sround(rad * points3)
    return(points3)
}
