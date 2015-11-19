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
  Xname <- formula[[2]]
  X <- eval(substitute(Xname), parent.frame())
  stopifnot(inherits(X, c("sp2", "sp3")))
  win <- X$win
  nX <- nrow(X$X)
  lam <- intensitysph(X)
  rho <- 2 * lam
  ## generate Poisson dummy points with same average intensity
  dum <- rpoispp.sp2(rho, win)
  df <- rbind(cbind(allcoords(X), isdata=TRUE),
              cbind(allcoords(dum), isdata=FALSE))
  df$rho <- rho
  fmla <- formula
  fmla[[2]] <- as.name("isdata")
  fit <- glm(fmla, data=df, family=binomial, offset=-log(rho))
  result <- list(rho=rho,
                 fit=fit)
  class(result) <- c("sphppm", class(result))
  return(result)
}

print.sphppm <- function(x, ...) {
  cat("Parametric model fitted to spherical point pattern data\n")
  fit <- x$fit
  cat("Model formula:", pasteFormula(formula(fit)), "\n")
  cat("Fitted coefficients:\n")
  print(coef(fit))
  return(invisible(NULL))
}

fitted.sphppm <- function(object, ...) {
  fit <- object$fit
  value <- with(fit$data,
                (rho * exp(predict(fit, type="link")))[isdata])
  return(value)
}

anova.sphppm <- function(object, ...) {
  anova(object$fit, ...)
}


