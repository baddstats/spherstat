Gstab.sphwin <- function (X, ...) {
    G <- Gsphere(X=X, ...)
    rad <- X$win$rad
    Gstab <- eval.fv(asin(sqrt(G))/(pi/2))
    Gstab <- rebadge.fv(Gstab, quote(Gstab(r)), "Gstab", names(G), new.labl = attr(G, 
        "labl"))
    return(Gstab)
}


Ginv.sphwin <- function(X, ...) {
   G <- Gsphere(X, ...)
   intX <- intensitysph(X)
   rad <- X$win$rad
   Ginv <- rad*acos(cround(sround(1 + (log(1+G)/(2*pi*intX*(rad^2))))))
   Ginv <- rebadge.fv(Ginv, quote(Ginv(r)), "Ginv", names(G), new.labl=attr(G,"labl"))
   return(Ginv)
}
