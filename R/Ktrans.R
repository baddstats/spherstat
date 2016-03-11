Kinv.sphwin <- function (X, ...) 
{
    K <- Ksphere(X, ...)
    rad <- X$win$rad
    Kinv <- eval.fv(acos(cround(sround(1-(K/(2*pi*(X$win$rad^2)))))))
    Kinv <- rebadge.fv(Kinv, quote(Kinv(r)), "Kinv", names(K), new.labl = attr(K, 
        "labl"))
    return(Kinv)
}

Kstab.sphwin <- function (X, ...) {
    K <- Ksphere(X, ...)
    Kstab <- eval.fv(sqrt(K))
    Kstab <- rebadge.fv(Kstab, quote(Kstab(r)), "Kstab", names(K), new.labl = attr(K, 
        "labl"))
    return(Kstab)
}
