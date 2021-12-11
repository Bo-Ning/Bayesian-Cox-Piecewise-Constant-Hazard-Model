#' Sampler for the intermediate intervals for the piecewise exponential prior
#' with dependent Gamma heights.
#'
#' This is the sampler for the intermediate intervals (= all intervals except
#' for the first and last one) in case the piecewise exponential prior with
#' dependent Gamma heights is selected. The sampler is described in the
#' Supplement to Castillo and Van der Pas (2020) and uses MCMC within Gibbs,
#' with a Gamma proposal with shape parameter equal to the number of events in
#' the interval plus some epsilon (to prevent proposals equal to zero if there
#' are no events in an interval) and rate parameter equal to the parameter alpha
#' (set by the user) divided by histogram height on the previous interval, plus
#' the total amount of time all individuals were exposed during this interval.
#' Most users of the package will not work with this function directly, but
#' instead use the main function \link{BayesSurv}, in which this particular
#' function is incorporated.
#'
#' @seealso \link{BayesSurv}, which computes the posterior mean and credible
#'   bands for the cumulative hazard and survival functions, as well as the
#'   posterior mean for the hazard. Within \link{BayesSurv}, the present
#'   function as well as \link{MCMCDepGammaFirst} is called through
#'   \link{SamplePosteriorDepGamma}.
#'
#' @references Castillo and Van der Pas (2020). Multiscale Bayesian survival
#'   analysis. <arXiv:2005.02889>.
#'
#' @param current The value of the height of the first interval from the
#'   previous iteration.
#' @param prev.haz The value of the height of the preceding interval from
#'   the previous iteration.
#' @param next.haz The value of the height of the next interval from the
#'   previous iteration.
#' @param failure The number of individuals who had an event during the first
#'   interval.
#' @param exposure The total amount of time all individuals were exposed for
#'   during the first interval.
#' @param alpha.dep The main parameter \eqn{\alpha} for the dependent Gamma
#'   prior, as described in the documentation for \link{BayesSurv}. It is
#'   recommended to take \code{alpha.dep} smaller than \code{alpha0.dep}.
#'
#' @return \item{res}{A new sample of the histogram height of the selected
#'   interval.}
#'
#' @export

MCMCDepGammaIntermediate<- function(current, prev.haz, next.haz, failure, exposure, alpha.dep){
  epsilon <- 0.1
  #The epsilon is to prevent proposals equal to zero (in case of
  #no failures in an interval)
  prop <- rgamma(1, failure + epsilon, alpha.dep/prev.haz + exposure)
  
  log.A <- alpha.dep*next.haz*(1/current - 1/prop) + epsilon*(log(current) - log(prop))

  if (is.na(log.A) == T) {
    log.A <- log(1e-6)
  }
  
  u <- runif(1)

  if(u < exp(log.A)){res <- prop}
  if(u >= exp(log.A)){res <- current}

  return(res)
}
