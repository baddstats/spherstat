#'
#' sphppm.R
#'
#'   Fit a parametric intensity model
#'   to point pattern data on a sphere
#'   and compute predicted intensity at each data point.
#'

allcoords <- function(X) {
  X2 <- if(inherits(X, "sp2")) X$X else convert2(X)
  X3 <- if(inherits(X, "sp2")) convert3(X) else X$X
  df <- cbind(X2, X3)
  if(!is.matrix(df)) df <- matrix(df, ncol=5)
  colnames(df) <- c("theta", "phi", "x1", "x2", "x3")
  return(as.data.frame(df))
}
    
sphppm <- function(formula) {
  callingframe <- parent.frame()
  Xname <- formula[[2]]
  X <- eval(substitute(Xname), parent.frame())
  stopifnot(inherits(X, c("sp2", "sp3")))
  win <- X$win
  nX <- nrow(X$X)
  lam <- intensitysph(X)
  rho <- 2 * lam
  ## generate Poisson dummy points with same average intensity
  dum <- rpoispp.sphwin(rho, win)
  df <- rbind(cbind(allcoords(X), isdata=TRUE),
              cbind(allcoords(dum), isdata=FALSE))
  df$rho <- rho
  ## change the LHS of the formula, and modify its environment
  fmla <- oldfmla <- formula
  fmla[[2]] <- as.name("isdata")
  fenv <- new.env(parent=callingframe)
  assign("df", df, envir=fenv)
  environment(fmla) <- fenv
  ## fit model
  fit <- glm(fmla, data=df, family=binomial, offset=-log(rho))
  ## pack up
  result <- list(rho=rho,
                 fit=fit,
                 X=X,
                 dum=dum,
                 oldfmla=oldfmla)
  class(result) <- c("sphppm", class(result))
  return(result)
}

print.sphppm <- function(x, ...) {
  cat("Parametric model fitted to spherical point pattern data\n")
  cat("Model formula:", pasteFormula(formula(x)), "\n")
  cat("Fitted coefficients:\n")
  print(coef(x))
  return(invisible(NULL))
}

coef.sphppm <- function(object, ...) {
  return(coef(object$fit))
}

formula.sphppm <- function(x, ...) { x$oldfmla }

fitted.sphppm <- function(object, ...) {
  fit <- object$fit
  value <- with(fit$data,
                (rho * exp(predict(fit, type="link")))[isdata])
  return(value)
}

anova.sphppm <- function(object, ...) {
  allargs <- list(object, ...)
  ismodel <- sapply(allargs, inherits, what="sphppm")
  allargs[ismodel] <- lapply(allargs[ismodel], getElement, name="fit")
  return(do.call(anova, allargs))
}

update.sphppm <- function(object, ...) {
  argh <- list(...)
  if(length(argh) == 0) return(object)
  if(any(ispat <- sapply(argh, inherits, what=c("sp2", "sp3")))) {
    ## One of the arguments is a point pattern 
    if(sum(ispat) > 1)
      stop("Only one point pattern should be given to update.sphppm",
           call.=FALSE)
    X <- argh[[ispat]]
    fenv <- environment(formula(object$fit))
    df <- get("df", fenv)
    dumdf <- df[with(df, !isdata), , drop=FALSE]
    df <- rbind(cbind(allcoords(X), isdata=TRUE, rho=object$rho),
                dumdf)
    assign("df", df, envir=fenv)
    assign("fmla", formula(object$fit), envir=fenv)
    newcall <- update(object$fit, evaluate=FALSE)
    newfit <- eval(substitute(newcall), envir=environment(formula(object$fit)))
    object$fit <- newfit
    object$X <- X
    if(all(ispat))
      return(object)
    ## additional arguments given: apply them now
    newobject <- do.call(update, append(list(object), argh[!ispat]))
    return(newobject)
  }
  newcall <- update(object$fit, ..., evaluate=FALSE)
  newfit <- eval(substitute(newcall), envir=environment(formula(object$fit)))
  object$fit <- newfit
  return(object)
}

predict.sphppm <- function(object, newdata=NULL,
                           type=c("intensity", "link", "terms"), ...) {
  fit <- object$fit
  type <- match.arg(type)
  if(is.null(newdata)) {
    value <- with(fit$data,
                  switch(type,
                         link = {
                           (log(rho) + predict(fit, type="link"))[isdata]
                         },
                         intensity = {
                           (rho * exp(predict(fit, type="link")))[isdata]
                         },
                         terms = {
                           tums <- predict(fit, type="terms")
                           con <- attr(tums, "constant")
                           val <- if(is.null(dim(tums))) tums[isdata] else 
                                  tums[isdata, , drop=FALSE]
                           attr(val, "constant") <- con
                           val
                         }))
    return(value)
  }
  if(!inherits(newdata, c("sp2", "sp3", "data.frame")))
    stop("newdata should be a data frame or a point pattern")
  if(inherits(newdata, c("sp2", "sp3"))) {
    newdata <- allcoords(newdata)
  } else {
    reqd <- c("theta", "phi", "x1", "x2", "x3") 
    if(!all(reqd %in% colnames(newdata)))
      stop(paste("newdata should have columns named", commasep(sQuote(reqd))))
  }
  rho1 <- with(fit$data, unique(rho))
  if(length(rho1) > 1)
    stop("Internal error: dummy intensity rho is not constant")
  newdata$rho <- rho1
  value <- with(fit$data,
                switch(type,
                       link = {
                         log(rho1) + predict(fit, newdata=newdata, type="link")
                         },
                         intensity = {
                           rho1 *
                             exp(predict(fit, newdata=newdata, type="link"))
                         },
                         terms = {
                           predict(fit, newdata=newdata, type="terms")
                         }))
  return(value)
}

simulate.sphppm <- function(object, nsim=1, ..., win, drop=TRUE) {
  if(missing(win) || is.null(win))
    win <- object$X$win
  lmax <- max(fitted(object),
              predict(object, newdata=object$dum, type="intensity"))
  obj <- object
  lambdafunction <- function(X, ..., fit=obj) {
    X <- sp2(X)
    predict(fit, newdata=X, type="intensity")
  }
  out <- replicate(nsim,
                   rpoispp.sphwin(lambda=lambdafunction,
                                  lmax=lmax,
                                  win = win),
                   simplify=FALSE)
  if(nsim == 1 && drop)
    out <- out[[1]]
  return(out)
}

is.poisson.sphppm <- function(x) TRUE

is.stationary.sphppm <- function(x) {
  trend <- rhs.of.formula(formula(x$fit))
  identical.formulae(trend, ~1)
}
