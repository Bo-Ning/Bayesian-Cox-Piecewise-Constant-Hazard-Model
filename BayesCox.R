#' Compute the posterior mean and a credible band for the survival and
#' cumulative hazard, and the posterior mean for the hazard.
#'
#' This is the main function of this package, computing relevant quantities for
#' a Bayesian survival analysis of (possibly right-censored) time-to-event-data.
#' Starting with a piecewise exponential prior with dependent or independent
#' Gamma heights (details below) on the hazard function, the function computes
#' the posterior mean for the hazard, cumulative hazard and survival function,
#' serving as an estimator for the true functions. In addition, for the
#' cumulative hazard and survival function, the radius for a fixed-width
#' credible band is computed. The interpretation of this credible band as a
#' confidence band is justified in Castillo and Van der Pas (2020).
#'
#' There are two options for the prior: a piecewise exponential (histogram)
#' prior with dependent Gamma heights and a piecewise exponential (histogram)
#' prior with independent Gamma heights. Both priors are described in detail in
#' Castillo and Van der Pas (2020). The dependent prior has a Markov structure,
#' where the height of each interval depends on the height of the previous
#' interval. It implements the autoregressive idea of Arjas and Gasbarra (1994).
#' With \eqn{\lambda_k} the histogram height on interval \eqn{k} and \eqn{\alpha} a
#' user-selected parameter, the structure is such that, with \eqn{K} the number
#' of intervals:
#'
#' \deqn{E[\lambda_k | \lambda_{k-1}, \ldots, \lambda_1] = \lambda_{k-1}, k = 2, \ldots, K.}
#' \deqn{Var(\lambda_k | \lambda_{k-1}, \ldots, \lambda_1) = (\lambda_{k-1})^2/\alpha, k = 2, \ldots, K.}
#'
#' In the independent Gamma prior, the prior draws for the \eqn{\lambda_k}'s are
#' independent of each other and are taken from a Gamma distribution with
#' user-specified shape and rate parameters.
#'
#' The guideline for the number of intervals \eqn{K} suggested by Castillo and
#' Van der Pas (2020) is
#'
#' \deqn{K = (n / \log{n})^{1/(1 + 2\gamma)},}
#'
#' where \eqn{n} is the number of observations and \eqn{\gamma} is related to
#' the smoothness of the true hazard function. In the absence of information
#' about the smoothness, a default value of \eqn{\gamma = 1/2} is recommended
#' and this is implemented as the default in this package. If this choice leads
#' to many intervals with zero events, it is recommended to decrease the number
#' of intervals.
#'
#' The samplers used for the dependent and independent Gamma priors are described
#' in the Supplement to Castillo and Van der Pas (2020).
#'
#'
#' @seealso \link{PlotBayesSurv} for a function that takes the result
#' from \code{BayesSurv()} and produces plots of the posterior mean of the hazard,
#' the posterior mean and credible band for the cumulative hazard, and the
#' posterior mean and credible band for the survival. To obtain direct samples
#' from the posterior for the hazard, see \link{SamplePosteriorDepGamma} and
#' \link{SamplePosteriorIndepGamma}.
#'
#' @references Castillo and Van der Pas (2020). Multiscale Bayesian survival
#'   analysis. <arXiv:2005.02889>.
#'
#' Arjas and Gasbarra (1994). Nonparametric Bayesian inference from right
#' censored survival data, using the Gibbs sampler. Statistica Sinica 4(2):505-524.
#'
#' @param df A dataframe, containing at minimum a column with follow-up times
#'   and a column with a status indicator (event observed or censored).
#' @param time The name of the column in the dataframe containing the (possibly
#'   right-censored) follow-up times, that is, the minimum of the time of the
#'   event and the time of censoring. Input the name as character/string.
#' @param event The name of the column in the dataframe containing the status
#'   indicator, which must be coded as: 0 = censored, 1 = event observed. Input
#'   the name as character/string.
#' @param prior Select either dependent or independent Gamma heights for the
#'   piecewise exponential prior on the hazard. Dependent heights (with the
#'   Markov structure described below) is default.
#' @param K The number of intervals to be used in the piecewise exponential
#'   (histogram) prior. Default is set to \eqn{K = (n / \log{n})^{1/2}}, with
#'   \eqn{n} the number of observations, as recommended by Castillo and Van der
#'   Pas (2020).
#' @param time.max The maximum follow-up time to consider, corresponding to the
#'   parameter \eqn{\tau} in Castillo and Van der Pas (2020).
#' @param alpha The function will compute (1-\code{alpha})100\% credible bands
#'   for the cumulative hazard and survival function.
#' @param N The number of samples to draw from the posterior.
#' @param alpha.dep For the dependent Gamma prior only. The main parameter
#'   \eqn{\alpha} for the dependent Gamma prior, as described below. It is
#'   recommended to take \code{alpha.dep} smaller than \code{alpha0.dep}.
#' @param alpha0.dep For the dependent Gamma prior only. The shape parameter for
#'   the Gamma prior on the histogram height for the first interval. It is
#'   recommended to take \code{alpha.dep} smaller than \code{alpha0.dep}.
#' @param beta0.dep For the dependent Gamma prior only. The rate parameter for
#'   the Gamma prior on the histogram height for the first interval.
#' @param alpha.indep For the independent Gamma prior only. The shape parameter
#'   for the Gamma prior on the histogram height for each interval.
#' @param beta.indep For the independent Gamma prior only. The rate parameter
#'   for the Gamma prior on the histogram height for each interval.
#' @param surv.factor The survival function is computed on an equispaced grid
#'  consisting of \code{K x surv.factor} (the number of intervals times this
#'   factor).
#' @param surv.epsilon The survival function is computed on the interval [0,
#'   \code{time.max} - \code{surv.epsilon}].
#'
#' @return \item{haz.post.mean}{The posterior mean for the hazard, given as the
#'   value on each of the \eqn{K} intervals.}
#'   \item{cumhaz.post.mean}{The posterior mean for the cumulative hazard,
#'   given as the value at the end of each of the \eqn{K} intervals. The
#'   cumulative hazard can be obtained from this by starting at 0 and
#'   linearly interpolating between each of the returned values.}
#'   \item{cumhaz.radius}{The radius for the credible set for
#'   the cumulative hazard.}
#'   \item{surv.post.mean}{The posterior mean for the
#'   survival, given at each value contained in the also returned
#'   \code{surv.eval.grid}.}
#'   \item{surv.radius}{The radius for the credible set for
#'   the survival.}
#'   \item{surv.eval.grid}{The grid on which the
#'   posterior mean for the survival has been computed. A finer grid can be
#'   obtained by increasing \code{surv.factor} in the function call.}
#'   \item{time.max}{The maximum follow-up time considered.}
#'
#' @examples #Demonstration on a simulated data set
#' library(simsurv)
#' library(ggplot2)
#' hazard.true <- function(t,x, betas, ...){1.2*(5*(t+0.05)^3 - 10*(t+0.05)^2 + 5*(t+0.05) ) + 0.7}
#' sim.df <- data.frame(id = 1:1000)
#' df <- simsurv(x = sim.df, maxt = 1, hazard = hazard.true)
#'
#' bs <- BayesSurv(df, "eventtime", "status")
#' PlotBayesSurv(bs, object = "survival")
#' PlotBayesSurv(bs, object = "cumhaz")
#' PlotBayesSurv(bs, object = "hazard")
#'
#' @export

BayesCox <- function(df, time = "time", event = "event",
                     prior = c("Dependent", "Independent"),
                     K = ceiling((dim(df)[1] / log(dim(df)[1]))^(1/2)),
                     time.max = max(df$time), alpha = 0.05,
                     N = 1000, alpha.dep = 1, alpha0.dep = 1.5, beta0.dep = 1,
                     alpha.indep = 1.5, beta.indep = 1, theta.prop.sd = NULL,
                     surv.factor = 10, surv.epsilon = 0.0000000001,
                     theta.int = NULL, lambda.int = NULL, z = NULL, 
                     burnin = NULL){
  
  # prior = match.arg(prior)
  
  # prior <- "Independent"
  # K <- ceiling((dim(df)[1] / log(dim(df)[1]))^(2/3))
  # time.max = max(df$time)
  # N <- 5000
  # alpha.indep = 1.5
  # beta.indep = 1
  
  if (is.vector(df$trt) == T) {
    p <- 1
  } else {
    p <- dim(df$trt)[2]
  }

  
  if (is.null(theta.int) == T) {
    model <- coxph(Surv(eventtime, status) ~ trt, data=df)
    theta.int <- model$coefficients # frequentist estimator
    if (p == 1) {
      exp.reg <- c(exp(df$trt * theta.int))
    } else {
      exp.reg <- c(exp(df$trt %*% theta.int))
    }
  } 
  
  if (is.null(theta.prop.sd) == T){
    theta.prop.sd <- sqrt(model$var)
  }
  
  if (is.null(lambda.int) == T) {
    bsline.cumhaz.coxph <- basehaz(model, centered = F) # true cumulative hazard
    interval <- c(0,cumsum(rep(1/K, K)))
    Lambda.int <- rep(NA, K)
    for (k in 1:K) {
      if (k == 1) {
        index <- bsline.cumhaz.coxph$time >= interval[k] & bsline.cumhaz.coxph$time <= interval[k+1]
        hzd.val <- bsline.cumhaz.coxph$hazard[index]
        time.val <- bsline.cumhaz.coxph$time[index]
        length.diff <- time.val - c(0, time.val[1:(sum(index)-1)])
        Lambda.int[k] <- sum(hzd.val * length.diff)/sum(length.diff)
      } else {
        index <- bsline.cumhaz.coxph$time >= interval[k] & bsline.cumhaz.coxph$time <= interval[k+1]
        if (sum(index) == 0) {
          Lambda.int[k] <- Lambda.int[k-1]
        } else {
          hzd.val <- bsline.cumhaz.coxph$hazard[index]
          time.val <- bsline.cumhaz.coxph$time[index]
          length.diff <- time.val - c(bsline.cumhaz.coxph$time[which(index == T)[1] - 1], time.val[1:(sum(index)-1)])
          Lambda.int[k] <- sum(hzd.val * length.diff)/sum(length.diff)
        }
      }
    }
    lambda.int <- Lambda.int - c(0,Lambda.int[1:(K-1)])
  }
  
  data.summary <- CoxReshapeData(df, time, event, exp.reg, lambda.int, K, time.max)
  fail <- data.summary$failures
  exp <- data.summary$exposures
  
  if(prior == "Dependent"){
    post.samples <- SamplePosteriorDepGamma(fail, exp, df, theta.int, lambda.int,
                                            p, N, alpha.dep, alpha0.dep, beta0.dep, theta.prop.sd)
  } #end Dependent

  if(prior == "Independent"){
    post.samples <- SamplePosteriorIndepGamma(fail, exp, df, theta.int, p,
                                              N, alpha.indep, beta.indep, theta.prop.sd)
  }
  
  if (is.null(burnin) == T) {
    burnin <- floor(N*0.2)
  }

  # collect result
  bsline.haz.samples <- post.samples$lambda.samples[(burnin+1):N, ]
  theta.samples <- post.samples$theta.samples[(burnin+1):N, ]
  if (p == 1) {
    theta.samples <- as.matrix(theta.samples)
  }
  
  # calculate quantiles
  # theta
  theta.pm <- apply(as.matrix(theta.samples), 2, mean)
  theta.lq <- quantile(theta.samples, probs = 0.025)
  theta.uq <- quantile(theta.samples, probs = 0.975)
  theta.radius <- RadiusCredibleSet(as.matrix(theta.samples), theta.pm, alpha)

  
  #baseline cumulative hazard function
  bs.haz.pm <- apply(bsline.haz.samples, 2, mean)
  bs.cumhaz <- t(apply(bsline.haz.samples*time.max/K, 1, cumsum))
  #divide by K to account for interval length
  bs.cumhaz.pm <- apply(bs.cumhaz, 2, mean)
  bs.cumhaz.radius <- RadiusCredibleSet(bs.cumhaz, bs.cumhaz.pm, alpha = alpha)
  # piecewise credible intervals
  bs.cumhaz.lq <- apply(bs.cumhaz, 2, quantile, probs = 0.025)
  bs.cumhaz.uq <- apply(bs.cumhaz, 2, quantile, probs = 0.975)
  bs.cumhaz.lb <- bs.cumhaz.pm - bs.cumhaz.radius
  bs.cumhaz.ub <- bs.cumhaz.pm + bs.cumhaz.radius

  
  if (is.null(z) == T) {z <- mean(df$trt)}
  #conditional cumulative hazard
  cond.cumhaz.samples <- 
    t(sapply(1:dim(bs.cumhaz)[1], 
            FUN = function(j) {bs.cumhaz[j, ] * exp(c(theta.samples[j,] %*% z))}))
  cond.cumhaz.pm <- apply(cond.cumhaz.samples, 2, mean)
  cond.cumhaz.radius <- RadiusCredibleSet(cond.cumhaz.samples, bs.cumhaz.pm, 
                               alpha = alpha)
  # piecewise credible intervals
  cond.cumhaz.lq <- apply(cond.cumhaz.samples, 2, quantile, probs = 0.025)
  cond.cumhaz.uq <- apply(cond.cumhaz.samples, 2, quantile, probs = 0.975)
  cond.cumhaz.lb <- cond.cumhaz.pm - cond.cumhaz.radius
  cond.cumhaz.ub <- cond.cumhaz.pm + cond.cumhaz.radius

  
  # baseline survival
  bs.survival <- t(apply(bs.cumhaz, 1,
                         SurvivalFromRegAndCumhaz, time.max = time.max,
                         surv.factor = surv.factor, surv.epsilon = surv.epsilon))
  bs.surv.pm <- apply(bs.survival, 2, mean)
  bs.surv.radius <- RadiusCredibleSet(bs.survival, bs.surv.pm, 
                                            alpha = alpha)
  # piecewise credible intervals
  bs.surv.lq <- apply(bs.survival, 2, quantile, probs = 0.025)
  bs.surv.uq <- apply(bs.survival, 2, quantile, probs = 0.975)
  # credible bands
  bs.surv.lb <- bs.surv.pm - bs.surv.radius
  bs.surv.ub <- bs.surv.pm + bs.surv.radius

  
  # conditional survival function
  
  #survival function, default is the baseline survival function
  # calculate the survivals, for any vector/scalar a, obtain the conditional survival function given z
  cond.surv <- t(apply(cond.cumhaz.samples, 1,
                       SurvivalFromRegAndCumhaz, time.max = time.max,
                       surv.factor = surv.factor, surv.epsilon = surv.epsilon))
  cond.surv.pm <- apply(cond.surv, 2, mean)
  cond.surv.radius <- RadiusCredibleSet(cond.surv, cond.surv.pm, 
                                            alpha = alpha)
  # piecewise credible intervals
  cond.surv.lq <- apply(cond.surv, 2, quantile, probs = 0.025)
  cond.surv.uq <- apply(cond.surv, 2, quantile, probs = 0.975)
  # credible bands
  cond.surv.lb <- cond.surv.pm - cond.surv.radius
  cond.surv.ub <- cond.surv.pm + cond.surv.radius
  
  
  # create frames
  cumhaz.eval.grid <- seq(0, time.max, length = dim(bs.cumhaz)[2])
  surv.eval.grid <- seq(0, time.max - surv.epsilon, length = dim(bs.survival)[2])
  theta.qtl <- data.frame(theta.pm, theta.lq, theta.uq, theta.radius)
  bs.cumhaz.qtl <- data.frame(bs.cumhaz.pm, bs.cumhaz.lq, bs.cumhaz.uq,
                              bs.cumhaz.lb, bs.cumhaz.ub, cumhaz.eval.grid)
  cond.cumhaz.qtl <- data.frame(cond.cumhaz.pm, cond.cumhaz.lq, cond.cumhaz.uq,
                                cond.cumhaz.lb, cond.cumhaz.ub, cumhaz.eval.grid)
  bs.surv.qtl <- data.frame(bs.surv.pm, bs.surv.lq, bs.surv.uq, bs.surv.lb, 
                            bs.surv.ub, surv.eval.grid)
  cond.surv.qtl <- data.frame(cond.surv.pm, cond.surv.lq, cond.surv.uq, 
                              cond.surv.lb, cond.surv.ub, surv.eval.grid)
  z.bar <- mean(df$trt)
  
  return(list(theta.qtl = theta.qtl, bs.cumhaz.qtl = bs.cumhaz.qtl,
              cond.cumhaz.qtl = cond.cumhaz.qtl, bs.surv.qtl = bs.surv.qtl,
              cond.surv.qtl = cond.surv.qtl, time.max = time.max, z.bar = z.bar)
  )
}
