Fstab.sphwin <- function (X, ...) {
    FF <- Fsphere(X, ...)
    rad <- X$win$rad
    Fstab <- eval.fv(asin(sround(cround(sqrt(FF))/(pi/2))))
    Fstab <- rebadge.fv(Fstab, quote(Fstab(r)), "Fstab", names(FF), new.labl = attr(FF, 
        "labl"))
    return(Fstab)
}

Finv.sphwin <- function(X, ...) {
   F <- Fsphere(X, ...)
   intX <- intensitysph(X)
   rad <- X$win$rad
   Finv <- rad*acos(cround(sround(1 + (log(1+F)/(2*pi*intX*(rad^2))))))
   Finv <- rebadge.fv(Finv, quote(Finv(r)), "Finv", names(F), new.labl=attr(F,"labl"))
   return(Finv)
}