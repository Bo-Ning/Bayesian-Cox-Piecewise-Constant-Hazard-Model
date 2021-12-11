#' Computes the radius of a fixed width credible set for the
#' survival or the cumulative hazard
#'
#' This function finds a radius such that (1-alpha)100\% of posterior
#' draws are within a distance of at most this radius to the posterior
#' mean. Most users will not use this function directly, but instead
#' use \link{BayesSurv}, in which this function is used.
#'
#' @seealso \link{BayesSurv}, which computes the posterior mean and the
#' radius of the credible band for the cumulative hazard function as
#' well as the survival, and the posterior mean for the hazard. These objects
#' can then be visualized by using \link{PlotBayesSurv}.
#'
#' @references Castillo and Van der Pas (2020). Multiscale Bayesian survival
#'   analysis. <arXiv:2005.02889>.
#'
#' @param draws A matrix of posterior draws of either the cumulative
#' hazard or the survival. Each row contains a draw, the columns correspond
#' to time points on which the cumulative hazard or survival is evaluated.
#'
#' @param post.mean The posterior mean of the cumulative hazard or survival
#' function, evaluated at the same time points as the draws.
#'
#' @param alpha The credible band will be such that (1-alpha)100\% of draws is
#' contained in it.
#'
#' @return \item{radius}{The radius of the credible set.}
#'
#' @export

RadiusCredibleSet <- function(draws, post.mean, alpha = 0.05){

  N <- dim(draws)[1]

  max.dist.per.draw <- apply(abs(sweep(draws, 2, post.mean)), 1, max)

  radius <- sort(max.dist.per.draw)[ceiling((1-alpha)*N)]

  return(radius)
}
