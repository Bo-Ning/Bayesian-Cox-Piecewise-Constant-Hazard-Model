#' Reshape right censored data to be used with a piecewise exponential
#' prior.
#'
#' To draw from the posterior of the piecewise exponential priors implemented
#' in this package, it is convenient to convert the data so that three vectors
#' are obtained: the first containing the total amount of time all individuals
#' were under follow-up during each interval, the second containing the
#' number of events that happened during each interval, and the third containing 
#' the number of events weighted by its exponential of regression coefficients 
#' in each interval. 
#' This function takes a dataframe with a column of times 
#' (the minimum of the time of the event
#' and the time of censoring) and a column indicating the status (0 if censored,
#' 1 if the event was observed) and reshapes it into the desired format. Most
#' users will not use this function directly, but will instead use the main
#' function \link{BayesSurv}, which uses the present function.
#'
#' @seealso \link{BayesSurv}, which computes the posterior mean and the
#' radius of the credible band for the cumulative hazard function as
#' well as the survival, and the posterior mean for the hazard. These objects
#' can then be visualized by using \link{PlotBayesSurv}.
#'
#' @references 1. Ning and Castillo (2021). Multiscale analysis of the Bayesian Cox 
#'   proportional hazard model <arXiv:2005.02889>.
#'
#' @param df A dataframe, containing at minimum a column with follow-up times
#'   and a column with a status indicator (event observed or censored).
#' @param time The name of the column in the dataframe containing the (possibly
#'   right-censored) follow-up times, that is, the minimum of the time of the
#'   event and the time of censoring. Input the name as character/string.
#' @param event The name of the column in the dataframe containing the status
#'   indicator, which must be coded as: 0 = censored, 1 = event observed. Input
#'   the name as character/string.
#' @param exp.reg The exponential of the regression term
#' @param K The number of intervals to be used in the piecewise exponential
#'   (histogram) prior. Default is set to \eqn{K = (n / \log{n})^{1/2}}, with
#'   \eqn{n} the number of observations, as recommended by Castillo and Van der
#'   Pas (2020).
#' @param time.max The maximum follow-up time to consider, corresponding to the
#'   parameter \eqn{tau} in Castillo and Van der Pas (2020).
#'
#' @return \item{failures}{A vector of length \eqn{K}, containing for each
#' interval the number of individuals who had an event during that interval.}
#' \item{exposures}{A vector of length \eqn{K}, containing for each interval
#' the total amount of time all individuals together were under follow-up
#' during that interval.}
#'
#' @export

CoxReshapeData <- function(df, time = "time", event = "event", exp.reg, lambda,
                           K = ceiling((dim(df)[1] / log(dim(df)[1]))^(1/2)),
                           time.max = max(df$time)){
  #Reshape the data
  #Input: dataframe with columns 'time' and 'event'
  #Output: number of failures and total exposure time per interval
  
  df.small <- data.frame(time = df$time, event = df$event, exp.reg = exp.reg)
  
  df.long <- survSplit( Surv(time, event) ~ exp.reg, data = df.small,
                        cut = seq(0, time.max, length.out = (K+1)),
                        id = "id", episode = "interval")
  
  df.long$interval <- (df.long$interval - 1) #interval count starts at 2 otherwise
  
  df.long$delta <- (df.long$time - df.long$tstart) * df.long$exp.reg  #time weighted by exp.reg spent per interval
  
  id.rep <- df.long$id
  
  # aggregate the failures, exposures, and exposures times e^(regression term)
  # according to each interval
  failures <- aggregate(event ~ interval, df.long, sum)$event
  exposures <- aggregate(delta ~ interval, df.long, sum)$delta
  
  #If nobody is under followup anymore during the final interval(s),
  #they will be discarded in the aggregate() - step. We augment the
  #data with zeroes for those intervals.
  indices.no.followup <- which( 1:K > max(df.long$interval) )
  failures[indices.no.followup] <- 0
  exposures[indices.no.followup] <- 0
  
  return(list(failures = failures, exposures = exposures))
}
