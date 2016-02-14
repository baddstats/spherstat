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

sphwin2owin <- function(w, eye, top) {
  ## project a spherical region to a 2D polygon
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
  type  <- w$type
  param <- w$param
  ref   <- w$ref
  curve1 <- curve2 <- NULL
  switch(type,
         sphere = { },
         band = ,
         bandcomp = {
           fullcircle <- seq(0, 2*pi, length=1000)
           if(type == "band") fullcircle <- rev(fullcircle)
           if(param[1] != 0)
             curve1 <- rot.sphere(cbind(param[1], fullcircle),
                                  northpole=ref, inverse=TRUE)
           if(param[2] != pi)
             curve2 <- rot.sphere(cbind(param[2], rev(fullcircle)),
                                  northpole=ref, inverse=TRUE)
         },
         wedge = {
           halfcircle <- seq(0, pi, length=500)
           long1 <- param[2]
           long2 <- long1 + param[1]
           curve1 <- rot.sphere(cbind(halfcircle, long1),
                                northpole=ref, inverse=TRUE)
           curve2 <- rot.sphere(cbind(rev(halfcircle), long2),
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
           colat <- sort(colat, decreasing=TRUE)
           long  <- sort(long)
           colats <- seq(colat[1], colat[2], length=250)
           longs  <- seq(long[1], long[2], length=250)
           curve1 <- rbind(cbind(colat[1], longs),
                           cbind(colats, long[2]),
                           cbind(colat[2], rev(longs)),
                           cbind(rev(colats), long[1]))
           curve1 <- rot.sphere(curve1, northpole=ref, inverse=TRUE)
         })
  V <- disc(3, npoly=8)
  if(!is.null(curve1)) {
    curve1 <- globelines(convert.globe(curve1), eye=eye, top=top, do.plot=FALSE)
    V <- intersect.owin(V, curve2bigowin(curve1))
  }
  if(!is.null(curve2)) {
    curve2 <- globelines(convert.globe(curve2), eye=eye, top=top, do.plot=FALSE)
    if(type == "bandcomp") {
      V <- union.owin(V, curve2bigowin(curve2))
    } else {
      V <- intersect.owin(V, curve2bigowin(curve2))
    }
  }
  intersect.owin(V, disc())
}

curve2bigowin <- function(a) {
  stopifnot(!is.null(dim(a)) && ncol(a) == 4)
  piece <- attr(a, "piece")
  Plist <- list()
  if(length(unique(piece)) > 1) {
    ## Although the original curve on the sphere was a single closed curve,
    ## it has been broken into several pieces. Reconnect them.
    aa <- split(as.data.frame(a), piece)
    starts <- lapply(aa, startof)
    ends <- lapply(aa, endof)
    edge.starts <- sapply(starts, on.edge)
    edge.ends <- sapply(ends, on.edge)
    closed <- edge.starts & edge.ends
    if(any(closed)) {
      ## curves which start and end at the edge of the world
      ## become separate polygons
      cc <- lapply(aa[closed], curve2bigowin)
      cc <- lapply(cc, getElement, name="bdry")
      Plist <- Reduce(append, cc)
      aa <- aa[!closed]
    }
    while(length(aa) > 1) {
      starts <- lapply(aa, startof)
      ends <- lapply(aa, endof)
      loose.starts <- !sapply(starts, on.edge)
      loose.ends <- !sapply(ends, on.edge)
      Starts <- as.ppp(do.call(concatxy,starts), W=square(c(-2,2)))
      Ends   <- as.ppp(do.call(concatxy,  ends), W=square(c(-2,2)))
      v <- nncross(Ends[loose.ends], Starts[loose.starts])
      hit <- (v$dist < 0.01)
      if(!any(hit))
        stop("Internal error: unable to repair broken pieces of curve")
      ipos <- which(hit)[1]
      jpos <- v$which[ipos]
      i <- which(loose.ends)[ipos]
      j <- which(loose.starts)[jpos]
      ## end of fragment i matches start of fragment j
      ## Coalesce them
      aa[[j]] <- rbind(aa[[i]], aa[[j]])
      aa <- aa[-i]
    }
    ## only one fragment remains
    a <- aa[[1]]
  }
  ## 'a' is a single unbroken curve
  n <- nrow(a)
  x <- c(a[,1], a[n,3])
  y <- c(a[,2], a[n,4])
  n <- n+1
  theta0 <- atan2(y[1], x[1])
  theta1 <- atan2(y[n], x[n])
  if(theta0 < 0) theta0 <- theta0 + 2 * pi
  if(theta1 < 1) theta1 <- theta1 + 2 * pi
  theta <- seq(theta1,
               if(theta1 < theta0) theta0 else (theta0+2*pi),
               length=64)
  P <- list(x=c(x, 3*cos(theta)),
            y=c(y, 3*sin(theta)))
  ## make polygonal owin
  Plist <- append(Plist, list(P))
  if(sum(sapply(Plist, Area.xypolygon)) < 0)
    P <- lapply(Plist, reverse.xypolygon)
  return(owin(poly=Plist))
}

on.edge <- function(xy, tol=0.001) (abs(sum(xy^2) - 1) < tol)
startof <- function(a) c(x=a[1,1], y=a[1, 2])
endof <- function(a)   c(x=a[nrow(a),3], y=a[nrow(a), 4])
