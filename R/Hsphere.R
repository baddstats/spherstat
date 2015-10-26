Hsphere <- function (X, refpoints, win=sphwin(type="sphere"), rvals=seq(from=0, to=pi, length=512), ...) {
    FF <- Fsphere(X, refpoints=refpoints, win=win, rvals=rvals, ...)
    rad <- X$win$rad
    Hsphere <- eval.fv(1-FF)
    if (any(varcols <- colnames(FF) %in% c("rip", "ls"))) {
        r <- with(Hsphere, .x)
        Hsphere[, varcols] <- as.data.frame(FF)[, varcols]/(2 * pi * (rad^2)* (1-cos(r/rad)))
        n <- npoints(X)
        A <- area.sphwin(X$win)
        if (any(colnames(FF) == "rip")) 
            Fstab[r == 0, "rip"] <- (2 * A/(n - 1)^2)/(4 * pi)
        if (any(colnames(FF) == "ls")) 
            Fstab[r == 0, "ls"] <- (2 * A/(n * (n - 1)))/(4 * pi)
    }
    Hsphere <- rebadge.fv(Hsphere, quote(Hsphere(r)), "Hsphere", names(FF), new.labl = attr(FF, 
        "labl"))
    return(Hsphere)
}
