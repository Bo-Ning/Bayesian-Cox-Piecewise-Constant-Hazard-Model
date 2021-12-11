#' Transform posterior draws from the cumulative hazard into posterior
#' draws from the survival.
#'
#' Most users will not use this function directly, but will instead use the main
#' function \link{BayesSurv}, which calls this function. This function may be
#' used to create posterior draws of the survival function, based on posterior
#' draws of the cumulative hazard. It does so at a number of equispaced time
#' points on the interval [0, \code{time.max} - \code{surv.epsilon}], with the
#' number equal to the product of the number of intervals used in the prior, and
#' a user-defined factor.
#'
#' @seealso \link{BayesSurv}, which computes the posterior mean and the
#' radius of the credible band for the cumulative hazard function as
#' well as the survival, and the posterior mean for the hazard. These objects
#' can then be visualized by using \link{PlotBayesSurv}.
#'
#' @references Castillo and Van der Pas (2020). Multiscale Bayesian survival
#'   analysis. <arXiv:2005.02889>.
#'
#' @param cumhaz A matrix containing posterior draws of the cumulative. Each row
#' contains one draw, the columns correspond to each interval. The values in each
#' draw are the values of the cumulative hazard at the end of the corresponding
#' interval.
#'
#' @param time.max The maximum follow-up time to consider, corresponding to the
#'   parameter \eqn{\tau} in Castillo and Van der Pas (2020).
#'
#' @param surv.factor The survival function is computed on an equispaced grid
#'   consisting of \code{K x surv.factor} (the number of intervals times this
#'   factor).
#'
#' @param surv.epsilon The survival function is computed on the interval [0,
#'   \code{time.max} - \code{surv.epsilon}].
#'
#' @return \item{surv(eval.vec)}{A numeric vector containing the posterior mean of
#' the survival function, evaluated at \code{K x surv.factor} (where \code{K}
#' is the number of intervals used in the prior) equidistant time points on
#' the interval [0, \code{time.max} - \code{surv.epsilon}].}
#'
#' @export

SurvivalFromRegAndCumhaz <- function(cumhaz.exp.reg, time.max, surv.factor = 10, surv.epsilon = 0.0000000001){
  #takes a vector of K values of the cumulative hazard at the intervals
  #evaluates the survival at a grid containing factor*K values
  #and returns this as a vector of length factor*K
  K <- length(cumhaz.exp.reg)
  eval.vec <- seq(0, (time.max - surv.epsilon), length = surv.factor*K)

  cumhaz.exp.reg.interpolated <- approxfun( (time.max/K)*(0:K), c(0, cumhaz.exp.reg) )
  surv <- function(t){exp(-cumhaz.exp.reg.interpolated(t))}

  return(surv(eval.vec))

}
