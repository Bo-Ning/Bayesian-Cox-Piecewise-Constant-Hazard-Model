#' Sampler for the first interval for the piecewise exponential prior
#' with dependent Gamma heights.
#'
#' This is the sampler for the first interval in case the piecewise exponential
#' prior with dependent Gamma heights is selected. The sampler is described in
#' the Supplement to Castillo and Van der Pas (2020). Most users of the package
#' will not work with this function directly, but instead use the main function
#' \link{BayesSurv}, in which this particular function is incorporated.
#'
#' @seealso \link{BayesSurv}, which computes the posterior mean and credible
#'   bands for the cumulative hazard and survival functions, as well as the
#'   posterior mean for the hazard. Within \link{BayesSurv}, the present
#'   function as well as \link{MCMCDepGammaIntermediate} is called through
#'   \link{SamplePosteriorDepGamma}.
#'
#' @references Castillo and Van der Pas (2020). Multiscale Bayesian survival
#'   analysis. <arXiv:2005.02889>.
#'
#' @param current The value of the height of the first interval from the
#'   previous iteration.
#' @param next.haz The value of the height of the second interval from the
#'   previous iteration.
#' @param failure The number of individuals who had an event during the first
#'   interval.
#' @param exposure The total amount of time all individuals were exposed for
#'   during the first interval.
#' @param alpha.dep The main parameter \eqn{\alpha} for the dependent Gamma
#'   prior, as described in the documentation for \link{BayesSurv}. It is
#'   recommended to take \code{alpha.dep} smaller than \code{alpha0.dep}.
#' @param alpha0.dep The shape parameter for the Gamma prior on the histogram
#'   height for the first interval. It is recommended to take \code{alpha.dep}
#'   smaller than \code{alpha0.dep}.
#' @param beta0.dep The rate parameter for the Gamma prior on the histogram
#'   height for the first interval.
#'
#' @return \item{res}{A new sample of the histogram height of the first interval.}
#'
#' @export

MCMCDepGammaFirst <- function(current, next.haz, failure, exposure, alpha.dep = 1, alpha0.dep = 1.5, beta0.dep = 1){
  #alpha0 > alpha.prior prevents proposals equal to zero
  prop <- rgamma(1, failure + alpha0.dep - alpha.dep, beta0.dep + exposure)

  log.A <- alpha.dep*next.haz*(1/current - 1/prop)
  u <- runif(1)

  if (is.na(log.A) == T) {
    log.A <- log(1e-6)
  }
  
  if(u < exp(log.A)){res <- prop}
  if(u >= exp(log.A)){res <- current}

  return(res)
}
