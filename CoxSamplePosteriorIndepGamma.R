#' Draw samples from the posterior for the regression and baseline hazard, using the piecewise
#' exponential (histogram) prior with independent Gamma heights
#'
#' The sampler is described in the Supplement to Ning and Castillo
#' (2021). Most users of the package will not work with this function directly,
#' but instead use the main function \link{BayesSurv}, in which this particular
#' function is incorporated.
#'
#' The samples returned by this function are draws from the posterior for the
#' hazard function. To obtain draws from the posterior for the cumulative
#' hazard, one can use numerical integration. One way to achieve this is by
#' first finding the values of the cumulative hazard at the end of each interval,
#' e.g. by \code{t(apply(samples*time.max/K, 1, cumsum))}, where \code{samples}
#' is the output from the present function and \code{time.max} and \code{K} are
#' as described for \link{BayesSurv}, and then using \code{approxfun()} to linearly
#' interpolate in between. To obtain posterior samples from the survival, one
#' could then use \link{SurvivalFromCumhaz}.
#'
#' @seealso \link{BayesSurv}, which computes the posterior mean and credible
#'   bands for the cumulative hazard and survival functions, as well as the
#'   posterior mean for the hazard. Within \link{BayesSurv}, the present
#'   function is called.
#'
#' @references 
#'
#' @param failures A vector of length \eqn{K} (the total number of intervals),
#'   containing for each interval the number of individuals who had an event
#'   during that interval.
#' @param exposures A vector of length \eqn{K} (the total number of intervals),
#'   containing for each interval the total amount of time all individuals
#'   together times their corresponding e^(theta'z) were under follow-up during that interval.
#' @param N The number of draws to take.
#' @param p The number of covariates (dimension of theta)
#' @param covs contains the value of covariates
#' @param theta.int initial value for theta for MCMC
#' @param alpha.indep The shape parameter for the Gamma prior on the histogram
#'   height for each interval.
#' @param beta.indep The rate parameter for the Gamma prior on the histogram
#'   height for each interval.
#' @param theta.prop.var The standard deviation for the proposal density of theta
#'   the proposal density is normally distributed centered at the previous draw 
#'   of theta
#' @return \item{samples}{A \eqn{N} by \eqn{K} (the number of draws by the
#' number of intervals) matrix, with each row containing a draw from the
#' posterior for the hazard, based on a histogram prior with independent
#' Gamma heights.}
#'
#' @export

SamplePosteriorIndepGamma <- 
  function(failures, exposures, df, theta.int, p = 1, N = 1000, 
           alpha.indep = 1.5, beta.indep = 1, theta.prop.var){
    
    K <- length(failures)
    theta.samples <- matrix(0, nrow = N, ncol = p)
    lambda.samples <- matrix(0, nrow = N, ncol = K)
    if (p > 1) {
      exp.reg <- c(exp(df$trt %*% theta.int))
    } else {
      exp.reg <- c(exp(df$trt * theta.int))
    }
    
    # start MCMC
    for (iter in 1:N) {  
      #### sample lambda #########
      lambda <- sapply(1:K, function(k) {
        rgamma(1, shape = (failures[k] + alpha.indep), rate = (exposures[k] + beta.indep))
      }) # sample baseline hazard function
      
      Lambda <- sapply(df$time, FUN = Lambda.at.time.t, lambda = lambda) 
      
      ##### sample theta #########
      # theta.prop <- rnorm(n = 1, mean = theta.int, sd = theta.prop.var)
      if (iter == 1) {
        theta <- theta.int
      }
      # theta.prop <- rnorm(n = 1, mean = theta, sd = theta.prop.var)
      if (p > 1) {
        theta.prop <- c(rmvnorm(n = 1, mean = theta, sigma = theta.prop.var))
        exp.reg.prop <- c(exp(df$trt %*% theta.prop))
      } else {
        theta.prop <- rnorm(n = 1, mean = theta, sd = sqrt(theta.prop.var))
        exp.reg.prop <- c(exp(df$trt * theta.prop))
      }
      
      # update baseline cumulative hazard function
      # data.summary.update <- CoxReshapeData(df = df, exp.reg = exp.reg.prop, lambda = lambda, K = K)
      # failures <- data.summary.update$failures
      # exposures <- data.summary.update$exposures
      
      
      # evaluate log posterior
      if (p > 1) {
        exp.reg <- c(exp(df$trt %*% theta))
      } else {
        exp.reg <- c(exp(df$trt * theta))
      }
      
      MH.logratio <- -sum(Lambda * (exp.reg.prop - exp.reg)) +
        sum(log((exp.reg.prop/exp.reg)[df$event == 1])) -
        sum(theta.prop^2/2) + sum(theta^2/2)
      # cat(MH.logratio, "\n")
      if(log(runif(1)) < MH.logratio){
        theta <- theta.prop       # accept move with probability min(1,A)
        exp.reg <- exp.reg.prop
      } else {
        if (iter == 1) {
          theta <- theta.int
        }
      }
      data.summary.update <- CoxReshapeData(df = df, exp.reg = exp.reg, lambda = lambda, K = K)
      failures <- data.summary.update$failures
      exposures <- data.summary.update$exposures
      
      #print(theta)
      theta.samples[iter, ] <- as.matrix(theta)
      lambda.samples[iter, ] <- lambda
    }
  
  return(list(lambda.samples = lambda.samples, theta.samples = theta.samples))
  
}
