\name{envelope.sphere}
\alias{envelope.sp2}
\alias{envelope.sphppm}
\title{Simulation Envelopes of Summary Function for Spherical Point Pattern}
\description{
  \code{envelope.sp2} computes simulation envelopes of a summary function
  for a point pattern on a sphere; \code{envelope.sphppm} computes simulation envelopes of a summary function,
  based on a fitted point process model, for a point pattern on a sphere.
}
\usage{
  \method{envelope}{sp2}(Y, fun=Ksphere, nsim=99, nrank=1, \dots, 
  funargs=list(), simulate=NULL, fix.n=FALSE, verbose=TRUE, 
  transform=NULL,global=FALSE,ginterval=NULL,use.theory=NULL,
  alternative=c("two.sided", "less", "greater"),
  scale=NULL, clamp=FALSE, 
  savefuns=FALSE, savepatterns=FALSE,
  nsim2=nsim, VARIANCE=FALSE, nSD=2, Yname=NULL, maxnerr=nsim,
  do.pwrong=FALSE, envir.simul=NULL)
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
    A spherical point pattern (object of class
    \code{"sp2"}).
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
    envelopes. Only relevant if \code{global=TRUE}.
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
\value{
  A function value table (object of class \code{"fv"})
  which can be plotted directly.
  See \code{\link[spatstat]{envelope}} for further details.
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
  The function \code{envelope.sp2} is the method for 
  spherical point patterns (objects of class \code{"sp2"}), and the function \code{envelope.sphppm} is the method for 
  point process models fitted to spherical point patterns
  (objects of class \code{"sphppm"}).

  For the most basic use, if you have a spherical point pattern \code{X} and
  you want to test Complete Spatial Randomness (CSR), type
  \code{plot(envelope(X, Ksphere,nsim=39))} to see the spherical
  \eqn{K} function for \code{X} plotted together with the envelopes of
  the spherical \eqn{K} function for 39 simulations of CSR. 

However, if you have fitted a model \code{M} to
  spherical point pattern data and you want to test the goodness-of-fit
  of the model, simply do \code{plot(envelope(M, nsim=39))}
  to display pointwise envelopes or
  \code{plot(envelope(M, nsim=19, global=TRUE))} for global envelopes.
  Both displays can be interpreted as performing a Monte Carlo test with
  significance level 0.05.
  
  To create simulation envelopes, the command \code{envelope(Y, ...)} 
  first generates \code{nsim} random point patterns
  in one of the following ways. 
  \itemize{
    \item 
    If \code{simulate=NULL},
    then we generate \code{nsim} simulations of
    Complete Spatial Randomness (i.e. \code{nsim} simulated point patterns
    each being a realisation of the uniform Poisson point process)
    with the same intensity as the pattern \code{Y}.
    \item
    If \code{simulate} is supplied, then it determines how the
    simulated point patterns are generated.
    See \code{\link{envelope}} for details.
  }
  
  The summary statistic \code{fun} is applied to each of these simulated
  patterns. Typically \code{fun} is one of the functions
  \code{\link{Ksphere}}, \code{\link{Gsphere}},
  \code{\link{Fsphere}} or \code{\link{Jsphere}}.
  It may also be a character string
  containing the name of one of these functions.
  
  For further information, see the documentation for
  \code{\link[spatstat]{envelope}}.
}
\references{
Baddeley, A., Diggle, P.J., Hardegen, A., Lawrence, T., Milne, R.K. and Nair, G. (2014) On tests of spatial pattern based on simulation envelopes. \emph{Ecological Monographs}, \bold{84}(3) 477--489

Cressie, N.A.C. \emph{Statistics for spatial data}. John Wiley and Sons, 1991. 

Diggle, P.J. \emph{Statistical analysis of spatial point patterns}. Arnold, 2003. 

Lawrence, T.J. (2017) Master's Thesis, University of Western Australia.

Ripley, B.D. (1981) \emph{Spatial statistics}. John Wiley and Sons. 

Ripley, B.D. \emph{Statistical inference for spatial processes}. Cambridge University Press, 1988. 

Stoyan, D. and Stoyan, H. (1994) \emph{Fractals, random shapes and point fields: methods of geometrical statistics}. John Wiley and Sons. 

} 
\seealso{

  \code{\link{sp2}},
  \code{\link{runif.sphwin}},
  \code{\link{rpoispp.sphwin}},
  \code{\link{Ksphere}},
  \code{\link{Fsphere}},
  \code{\link{Gsphere}},
  \code{\link{Jsphere}}.
}
\examples{
   X <- runif.sphwin(30)
   # envelopes of Ripley's K for this pattern
   EX <- envelope(X, Ksphere, nsim=39)
   plot(EX)
   M <- sphppm(X ~ x1+x2+x3)
   # envelopes of Ripley's K for this model
   EM <- envelope(M, Ksphere, nsim=19)
   plot(EM) 
   # envelopes of inhomogeneous K, using correctly-updated intensities
   EM2 <- envelope(M, Ksphere, lambdavalues=M, nsim=19)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{htest}
\keyword{hplot}
\keyword{iteration}
 
 
