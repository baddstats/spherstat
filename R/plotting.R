## plot methods for sp2, sp3, sphwin

plot.sp3 <- plot.sp2 <- function(x, ..., eye, top, add=FALSE) {
  X <- convert.globe(x)
  if(missing(eye) || is.null(eye)) {
    eye <- place('nedlands')
  } else if(!is.globe.point(eye)) {
    eye <- Convert.globe(eye)
  }
  if(missing(top) || is.null(top)) {
    top <- place('northpole')
  } else if(!is.globe.point(top)) {
    top <- Convert.globe(top)
  }
  if(!add) {
    globeearth(NULL)
    plot(x$win, ..., eye=eye, top=top, add=TRUE)
  }
  globepoints(X, eye=eye, top=top, ...)
}

plot.sphwin <- function(x, ..., eye, top, add=FALSE) {
  if(!add)
    globeearth(NULL)
  if(missing(eye) || is.null(eye)) {
    eye <- place('nedlands')
  } else if(!is.globe.point(eye)) {
    eye <- Convert.globe(eye)
  }
  if(missing(top) || is.null(top)) {
    top <- place('northpole')
  } else if(!is.globe.point(top)) {
    top <- Convert.globe(top)
  }
  type  <- x$type
  param <- x$param
  ref   <- x$ref
  curve1 <- curve2 <- NULL
  switch(type,
         sphere = { },
         band = ,
         bandcomp = {
           fullcircle <- seq(0, 2*pi, length=1000)
           if(param[1] != 0)
             curve1 <- rot.sphere(cbind(param[1], fullcircle),
                                  northpole=ref, inverse=TRUE)
           if(param[2] != pi)
             curve2 <- rot.sphere(cbind(param[2], fullcircle),
                                  northpole=ref, inverse=TRUE)
         },
         wedge = {
           halfcircle <- seq(0, pi, length=500)
           long1 <- param[2]
           long2 <- long1 + param[1]
           curve1 <- rot.sphere(cbind(halfcircle, long1),
                                northpole=ref, inverse=TRUE)
           curve2 <- rot.sphere(cbind(halfcircle, long2),
                                northpole=ref, inverse=TRUE)
         },
         polygon = {
           verts <- param
           curve1 <- matrix(, 0, 2)
           nv <- nrow(verts) - 1
           for(i in 1:nv) {
             path <- geodesicarc(verts[i,], verts[i+1,])
             curve1 <- rbind(curve1, path)
           }
         },
         quadrangle = {
           colat <- param[1:2]
           long <- param[4] + c(0, param[3])
           colats <- seq(colat[1], colat[2], length=250)
           longs  <- seq(long[1], long[2], length=250)
           curve1 <- rbind(cbind(colat[1], longs),
                           cbind(colats, long[2]),
                           cbind(colat[2], rev(longs)),
                           cbind(rev(colats), long[1]))
           curve1 <- rot.sphere(curve1, northpole=ref, inverse=TRUE)
         })
  if(!is.null(curve1))
    globelines(convert.globe(curve1), ..., eye=eye, top=top)
  if(!is.null(curve2))
    globelines(convert.globe(curve2), ..., eye=eye, top=top)
  return(invisible(NULL))
}

is.globe.point <- function(x) {
  is.list(x) && identical(names(x), c("lon", "lat"))
}

Convert.globe <- function(x) {
  if(!is.matrix(x)) x <- matrix(x, nrow=1)
  convert.globe(x)
}

geodesicarc <- function(A, B, n=100) {
  D <- rot.sphere(A, northpole=B)
  colats <- seq(D[1], 0, length=n)
  Path <- cbind(colats, D[2])
  P <- rot.sphere(Path, northpole=B, inverse=TRUE)
  return(P)
}
