# methods for 'rotate' for spherical objects

rotate.sphwin <- function(X, northpole, inverse=FALSE, ...) {
  if(length(list(...)) > 0)
    warning("Additional arguments were ignored by rotate.sphwin")
  switch(X$type,
         sphere={ },
         band = ,
         bandcomp = ,
         wedge = ,
         quadrangle = {
           ## rotate reference axis
           X$ref <- rot.sphere(X$ref,
                               northpole=northpole, inverse=inverse,
                               rad=X$rad)
         },
         polygon = {
           ## rotate polygon vertices
           X$param <- rot.sphere(X$param,
                               northpole=northpole, inverse=inverse,
                               rad=X$rad)
         },
         stop("Unsupported window type"))
  return(X)
}

rotate.sp2 <- function(X, northpole, inverse=FALSE, ...) {
  if(length(list(...)) > 0)
    warning("Additional arguments were ignored by rotate.sp2")
  sp2(rot.sphere(X, northpole=northpole, inverse=inverse),
      rotate(X$win, northpole=northpole, inverse=inverse))
}

rotate.sp3 <- function(X, northpole, inverse=FALSE, ...) {
  if(length(...) > 0)
    warning("Additional arguments were ignored by rotate.sp3")
  sp3(rot.sphere(X, northpole=northpole, inverse=inverse),
      rotate(X$win, northpole=northpole, inverse=inverse))
}
