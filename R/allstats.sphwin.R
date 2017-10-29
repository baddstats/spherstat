allstats.sphwin <- function (X, ..., dataname = NULL, verbose = FALSE) 
{
    stopifnot(inherits(X, c("sp2", "sp3", "matrix"))) 
    if (verbose) 
        cat("Calculating F, G, J ...")
    Jout <- do.call.matched(Jsphere, list(X = X, ...))
    if (verbose) 
        cat("ok.\n")
    Fout <- attr(Jout, "F")
    Gout <- attr(Jout, "G")
    attr(Jout, "F") <- NULL
    attr(Jout, "G") <- NULL
    fns <- list(`F function` = Fout, `G function` = Gout, `J function` = Jout)
    if (verbose) 
        cat("Calculating K function...")
    Kout <- do.call.matched(Ksphere, list(X = X, ...))
    fns <- append(fns, list(`K function` = Kout))
    if (verbose) 
        cat("done.\n")
    if (is.null(dataname)) 
        dataname <- short.deparse(substitute(X))
    title <- paste("Four summary functions for ", dataname, ".", 
        sep = "")
    attr(fns, "title") <- title
    fns <- as.anylist(fns)
    return(fns)
}
