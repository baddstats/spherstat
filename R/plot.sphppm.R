#' plot an intensity surface on the sphere

plot.sphppm <- function(x, eye=place("nedlands"), top=place("northpole"),
                        w,
                       ..., eps=NULL, dimyx=NULL, main="", 
                       action=c("image", "contour", "imagecontour"),
                       col.image=NULL, col.lines=NULL) {
  action <- match.arg(action)
  if(inherits(x, "sphppm")) {
    #' make grid of points for function evaluation
    lonlat <- expand.grid(lon=0:359,
                          lat=asin(seq(-0.999, 0.999, length=200)) * 180/pi)
    lonlat <- rbind(lonlat,
                    data.frame(lon=0, lat=c(90, -90)))
    lon <- lonlat$lon
    lat <- lonlat$lat
    phi <- (lon/180)*pi
    theta <- ((90-lat)/180)*pi
    Xgrid <- sp2(cbind(theta=theta,phi=phi))
    #' restrict to window
    if(missing(w)) w <- x$X$win
    if(!is.null(w)) {
      inside <- in.W(Xgrid, w)
      Xgrid$X <- Xgrid$X[inside,,drop=FALSE]
      Xgrid$win <- w
    }
    #' compute fitted intensity
    values <- predict(x, newdata=Xgrid)
    df <- cbind(Convert.globe(Xgrid$X),
                data.frame(values=values))
  } else if(is.data.frame(x)) {
    stopifnot(ncol(x) == 3)
    if(is.null(cnames <- colnames(x)) ||
       !identical(cnames[1:2], c("lon", "lat")))
      warning(paste("Interpreting the first two columns of data as",
                    "longitude and latitude respectively, in degrees"),
              call.=FALSE)
    df <- x
    if(missing(w)) w <- NULL
  } else stop("x should be a fitted model (sphppm)")

  if(!is.globe.point(eye)) eye <- Convert.globe(eye)
  if(!is.globe.point(top)) eye <- Convert.globe(top)
  eye3 <- ensure3d(eye)
  top3 <- ensure3d(top)
  
  spos <- spatialpos(df[,1], df[,2])
  mpos <- orthogproj(eye3, top3, spos)
  xx <- mpos[,1]
  yy <- mpos[,2]
  ok <- (mpos[,3] < 0)

  D <- disc(1)
  W <- disc(1, mask=TRUE, eps=eps, dimyx=dimyx)

  if(!is.null(w)) {
    wow <- sphwin2owin(w, eye=eye, top=top)
    W <- intersect.owin(W, wow)
  }
  X <- ppp(xx[ok], yy[ok], marks=df[ok,3], window=W, check=FALSE)
  ##  sigma <- 2 * mean(nndist(X))
  sigma <- 0.0125
  Y <- Smooth(X, sigma=sigma, eps=eps, dimyx=dimyx)

  switch(action,
         image = {
           do.call(plot.im,
                   resolve.defaults(list(x=Y, col=col.image),
                                    list(...),
                                    list(main=main, box=FALSE),
                                    .MatchNull=FALSE, .StripNull=TRUE))
           plot(D, add=TRUE)
         },
         contour = {
           plot(D, main=main)
           do.call(contour.im,
                   resolve.defaults(list(x=Y, col=col.lines, add=TRUE),
                                    list(...),
                                    .MatchNull=FALSE, .StripNull=TRUE))
         },
         imagecontour = {
           do.call.matched(plot.im,
                           resolve.defaults(list(x=Y, col=col.image),
                                            list(...),
                                            list(main=main, box=FALSE),
                                            .MatchNull=FALSE, .StripNull=TRUE),
                           extrargs="box")
           extrargs <- setdiff(names(formals(contour.default)),
                               union(c("x", "y", "z", "..."),
                                     names(formals(contour.im))))
           do.call.matched(contour.im,
                           resolve.defaults(list(x=Y, col=col.lines, add=TRUE),
                                            list(...),
                                            .MatchNull=FALSE, .StripNull=TRUE),
                           extrargs=extrargs)
           plot(D, add=TRUE)
         })
  if(!is.null(w)) plot(wow, add=TRUE)
  return(invisible(NULL))
}

