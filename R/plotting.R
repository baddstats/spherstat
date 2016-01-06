## plot methods for sp2, sp3, sphwin

plot.sp3 <- plot.sp2 <- function(x, ..., eye, top, add=FALSE) {
  X <- convert.globe(x)
  if(!add) {
    globeearth(NULL)
    plot(x$win, ..., add=TRUE)
  }
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
                                  northpole=ref)
           if(param[2] != pi)
             curve2 <- rot.sphere(cbind(param[2], fullcircle),
                                  northpole=ref)
         },
         wedge = ,
         polygon = ,
         quadrangle = {
           warning(paste(
             "Plotting is not yet implemented for windows of type",
             sQuote(type)))
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



