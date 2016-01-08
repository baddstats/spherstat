\name{envelope.sphppm}
\alias{envelope.sphppm}
\title{Simulation Envelopes for Spherical Point Process Model}
\description{
  Computes simulation envelopes of a summary function,
  based on a fitted point process model, for a point pattern on a sphere.
}
\usage{
\method{envelope}{sphppm}(Y, fun = Ksphere, nsim = 99, nrank = 1, \dots,
                 funargs = list(), simulate = NULL, fix.n = FALSE,
                 verbose = TRUE, transform = NULL, global = FALSE,
                 ginterval = NULL, use.theory = NULL,
                 alternative = c("two.sided", "less", "greater"),
                 scale = NULL, clamp = FALSE,
                 savefuns = FALSE, savepatterns = FALSE,
                 nsim2 = nsim, VARIANCE = FALSE, nSD = 2,
                 Yname = NULL, maxnerr = nsim, do.pwrong = FALSE,
                 envir.simul = NULL)
}
\arguments{
  \item{Y}{
    Fitted point process model on a sphere (object of class
    \code{"sphppm"} usually obtained from \code{\link{sphppm}}).
  }
  \item{fun}{
    Function that computes the desired summary statistic
    for a spherical point pattern. 
  }
  \item{nsim}{
    Number of simulated point patterns to be generated
    when computing the envelopes.
  }
  \item{nrank}{
    Integer. Rank of the envelope value amongst the \code{nsim} simulated
    values. A rank of 1 means that the minimum and maximum
    simulated values will be used.
  }
  \item{\dots}{
    Extra arguments passed to \code{fun}.
  }
  \item{funargs}{
    A list, containing extra arguments to be passed to \code{fun}.
  }
  \item{simulate}{
    Optional. Specifies how to generate the simulated point patterns.
    If \code{simulate} is an expression in the R language, then this
    expression will be evaluated \code{nsim} times,
    to obtain \code{nsim} point patterns which are taken as the
    simulated patterns from which the envelopes are computed.
    If \code{simulate} is a list of point patterns, then the entries
    in this list will be treated as the simulated patterns from which
    the envelopes are computed.
    Alternatively \code{simulate} may be an object produced by the
    \code{envelope} command: see Details.
  }
  \item{fix.n}{
    Logical flag indicating whether to generate simulations with the
    same number of points as in \code{Y}.
    (Ignored if the argument \code{simulate} is given.)
    (CURRENTLY IGNORED IN ALL CASES)
  }
  \item{verbose}{
    Logical flag indicating whether to print progress reports
    during the simulations.
  }
  \item{transform}{
    Optional. A transformation to be applied to the
    function values, before the envelopes are computed.
    An expression object (see Details).
  }
  \item{global}{
    Logical flag indicating whether envelopes should be pointwise
    (\code{global=FALSE}) or simultaneous (\code{global=TRUE}).
  }
  \item{ginterval}{
    Optional.
    A vector of length 2 specifying
    the interval of \eqn{r} values for the simultaneous critical
    envelopes. Only relevant if \code{global=TRUE}.
  }
  \item{use.theory}{
    Logical value indicating whether to use the theoretical value,
    computed by \code{fun}, as the reference value for simultaneous
    envelopes. Applicable only when \code{global=TRUE}.
  }
  \item{alternative}{
    Character string determining whether the envelope corresponds
    to a two-sided test (\code{side="two.sided"}, the default)
    or a one-sided test with a lower critical boundary
    (\code{side="less"}) or a one-sided test
    with an upper critical boundary (\code{side="greater"}).
  }
  \item{scale}{
    Optional. Scaling function for global envelopes.
    A function in the \R language which determines the
    relative scale of deviations, as a function of
    distance \eqn{r}, when computing the global envelopes.
    Applicable only when \code{global=TRUE}.
    Summary function values for distance \code{r}
    will be \emph{divided} by \code{scale(r)} before the
    maximum deviation is computed. The resulting global envelopes
    will have width proportional to \code{scale(r)}. 
  }
  \item{clamp}{
    Logical value indicating how to compute envelopes when
    \code{alternative="less"} or \code{alternative="greater"}.
    Deviations of the observed
    summary function from the theoretical summary function are initially
    evaluated as signed real numbers, with large positive values indicating
    consistency with the alternative hypothesis.
    If \code{clamp=FALSE} (the default), these values are not changed.
    If \code{clamp=TRUE}, any negative values are replaced by zero.
  }
  \item{savefuns}{
    Logical flag indicating whether to save all the simulated
    function values.
  }
  \item{savepatterns}{
    Logical flag indicating whether to save all the simulated
    point patterns.
  }
  \item{nsim2}{
    Number of extra simulated point patterns to be generated
    if it is necessary to use simulation to estimate the theoretical
    mean of the summary function. Only relevant when \code{global=TRUE}
    and the simulations are not based on CSR.
  }
  \item{VARIANCE}{
    Logical. If \code{TRUE}, critical envelopes will be calculated
    as sample mean plus or minus \code{nSD} times sample standard
    deviation.
  }
  \item{nSD}{
    Number of estimated standard deviations used to determine
    the critical envelopes, if \code{VARIANCE=TRUE}.
  }
  \item{Yname}{
    Character string that should be used as the name of the 
    data point pattern \code{Y} when printing or plotting the results.
  }
  \item{maxnerr}{
    Maximum number of rejected patterns.
    If \code{fun} yields an error when applied to a simulated point
    pattern (for example, because the pattern is empty and \code{fun}
    requires at least one point), the pattern will be rejected
    and a new random point pattern will be generated. If this happens
    more than \code{maxnerr} times, the algorithm will give up.
  }
  \item{do.pwrong}{
    Logical. If \code{TRUE}, the algorithm will also estimate
    the true significance level of the \dQuote{wrong} test (the test that
    declares the summary function for the data to be significant
    if it lies outside the \emph{pointwise} critical boundary at any
    point). This estimate is printed when the result is printed.
  }
  \item{envir.simul}{
    Environment in which to evaluate the expression \code{simulate},
    if not the current environment.
  }
}
\details{
  The \code{\link[spatstat]{envelope}} command performs simulations and
  computes envelopes of a summary statistic based on the simulations.
  The result is an object that can be plotted to display the envelopes.
  The envelopes can be used to assess the goodness-of-fit of
  a point process model to point pattern data.
  
  The \code{\link[spatstat]{envelope}} function is generic, with methods for
  the classes \code{"ppp"}, \code{"ppm"} and \code{"kppm"}
  described in the help file for \code{\link[spatstat]{envelope}}.
  This function \code{envelope.sphppm} is the method for 
  point process models fitted to spherical point patterns
  (objects of class \code{"sphppm"}).

  For the most basic use, if you have fitted a model \code{M} to
  spherical point pattern data and you want to test the goodness-of-fit
  of the model, simply do \code{plot(envelope(M, nsim=39))}
  to display pointwise envelopes or
  \code{plot(envelope(M, nsim=19, global=TRUE))} for global envelopes.
  Both displays can be interpreted as performing a Monte Carlo test with
  significance level 0.05.
}
\value{
  A function value table (object of class \code{"fv"})
  which can be plotted directly.
  See \code{\link[spatstat]{envelope}} for further details.
}
\references{
  Lawrence, T.J. (2016) Master's Thesis, University of Western Australia.
}
\author{
Adrian Baddeley.
}
\seealso{
  \code{\link{envelope.sp2}} for tests of complete spatial randomness.

  \code{\link{rpoispp.sphwin}},
  \code{\link{Ksphere}},
  \code{\link{Fsphere}},
  \code{\link{Gsphere}},
  \code{\link{Jsphere}}.
}
\examples{
   X <- runif.sphwin(30)
   M <- sphppm(X ~ x1+x2+x3)
   # envelopes of Ripley's K
   E <- envelope(M, Ksphere, nsim=19)
   plot(E) 
   # envelopes of inhomogeneous K, using correctly-updated intensities
   E2 <- envelope(M, Ksphere, lambdavalues=M, nsim=19)
}
\keyword{spatial}
\keyword{htest}
\keyword{hplot}
\keyword{iteration}
