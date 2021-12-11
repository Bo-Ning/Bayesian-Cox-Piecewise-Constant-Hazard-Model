#' Plot the posterior mean with credible band for the survival function
#' or cumulative hazard, or the posterior mean for the hazard
#'
#' This function takes the output from \link{BayesSurv} and uses ggplot2
#' to make plots of (1) the posterior mean of the survival function
#' with credible band, or (2) the posterior mean of the cumulative hazard
#' with credible band, or (3) the posterior mean of the cumulative hazard.
#' Users can select some plotting options within this function. Further
#' changes to the plot can be made by storing the plot and adding
#' ggplot2 syntax (see the examples).
#'
#' The posterior mean of the hazard and the posterior mean and credible
#' band of the cumulative hazard are plotted exactly. The survival is plotted
#' exactly at the points contained in the vector \code{surv.eval.grid} contained
#' in the object created by \link{BayesSurv}. Between these points, the survival
#' is linearly interpolated. To evaluate the survival exactly at more points
#' (and to obtain a smoother plot), increase the parameter \code{surv.factor}
#' within \link{BayesSurv}.
#'
#' @seealso \link{BayesSurv} to create the required object for this
#' plotting function.
#'
#' @references Castillo and Van der Pas (2020). Multiscale Bayesian survival
#'   analysis. <arXiv:2005.02889>.
#'
#' @param bayes.surv.object The output from the function \link{BayesSurv}.
#' @param object The object to be plotted, the user may select "survival"
#' for the survival function, "cumhaz" for the cumulative hazard, or
#' "hazard" for the hazard function. Default is the survival function.
#' @param band Indicator whether a credible band should be plotted (only
#' possible for the survival function and the cumulative hazard).
#' @param color The color to be used for the posterior mean and the
#' credible band (if applicable).
#' @param plot.title A title for the plot.
#' @param xlab A label for the horizontal axis.
#' @param ylab a label for the vertical axis.
#' @param legend If TRUE, a legend saying 'Credible band' will be included.
#' @param alpha.band The transparency of the credible band.
#'
#' @return \item{gg}{The plot, which may be edited further by adding
#' ggplot2 syntax.}
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
#'
#' cumhaz.plot <- PlotBayesSurv(bs, object = "cumhaz")
#' cumhaz.plot + labs(title = "Cumulative hazard")
#'
#' @export

PlotBayesSurv <- function(bayes.surv.object,
                            object = c("survival", "cumhaz", "hazard"),
                            band = TRUE, color = "darkblue",
                            plot.title = "", xlab = "time", ylab = "",
                            legend = TRUE, alpha.band = 0.4
                            ){

  object = match.arg(object)

  bs <- bayes.surv.object

  if(object == "survival"){
    surv.pm <- approxfun(bs$surv.eval.grid, bs$surv.post.mean)

    surv.pm.at.grid <- surv.pm(bs$surv.eval.grid)

    surv.radius <- bs$surv.radius

    df.surv <- data.frame(t = bs$surv.eval.grid,
                          post.mean = surv.pm.at.grid,
                          lower = pmax(0, surv.pm.at.grid - surv.radius),
                          upper = pmin(1, surv.pm.at.grid + surv.radius))

    gg <- ggplot(data = df.surv, aes(x = t, y = df.surv$post.mean, group = 1))
    gg <- gg + geom_line(colour = color, size = 1) + theme_classic()
    if(band == T){ gg <- gg + geom_ribbon(aes(ymin = df.surv$lower, ymax = df.surv$upper,
                               fill = "Credible band" ), alpha = alpha.band) }
    gg <- gg + theme(plot.title = element_text(size = 16))
    gg <- gg + theme(axis.text.x = element_text(size = 12))
    gg <- gg + theme(axis.text.y = element_text(size = 12))
    gg <- gg + scale_fill_manual("",values= c("Credible band" = color))
    gg <- gg + theme(legend.position = c(.75, 1),
                     legend.text = element_text(size = 12))
    if(legend == F){gg <- gg + theme(legend.position = "none")}
    gg <- gg + ylim(c(0, 1))
    gg <- gg + labs(title = plot.title)  + xlab(xlab) + ylab(ylab)
    #print(gg)
    return(gg)
  } #end survival

  if(object == "cumhaz"){
    K <- length(bs$haz.post.mean)

    cumhaz.pm <- approxfun(c(0, (bs$time.max/K)*(1:K) ),
                           c(0, cumsum(bs$haz.post.mean*bs$time.max/K)))

    cumhaz.pm.at.grid <- cumhaz.pm(bs$surv.eval.grid)

    cumhaz.radius <- bs$cumhaz.radius


    df.cumhaz <- data.frame(t = bs$surv.eval.grid,
                            post.mean = cumhaz.pm.at.grid,
                            lower = pmax(0, cumhaz.pm.at.grid - cumhaz.radius),
                            upper = cumhaz.pm.at.grid + cumhaz.radius)

    gg <- ggplot(data = df.cumhaz, aes(x = t, y = df.cumhaz$post.mean, group = 1))
    gg <- gg + geom_line(colour = color, size = 1) + theme_classic()
    if(band == T){gg <- gg + geom_ribbon(aes(ymin = df.cumhaz$lower, ymax = df.cumhaz$upper,
                               fill = "Credible band" ), alpha = alpha.band)}
    gg <- gg + theme(plot.title = element_text(size = 16))
    gg <- gg +   theme(plot.subtitle = element_text(size = 16))
    gg <- gg + theme(axis.text.x = element_text(size = 12))
    gg <- gg + theme(axis.text.y = element_text(size = 12))
    gg <- gg + scale_fill_manual("",values= c("Credible band" = color))
    gg <- gg + theme(legend.position = c(.2, 0.9),
                     legend.text = element_text(size = 12))
    if(legend == F){gg <- gg + theme(legend.position = "none")}
    gg <- gg + labs(title = plot.title)  + xlab(xlab) + ylab(ylab)
    #print(gg)
    return(gg)
  } #end cumhaz

  if(object == "hazard"){
    K <- length(bs$haz.post.mean)

    hazard.pm <- stepfun(seq(0, bs$time.max, length.out = (K+1)), c(0, bs$haz.post.mean, 0))

    df.hazard <- data.frame(t = bs$surv.eval.grid,
                            hazard = hazard.pm(bs$surv.eval.grid)
    )

    gg <- ggplot(data = df.hazard, aes(x = t, y = df.hazard$hazard, group = 1))
    gg <- gg + geom_point(colour = color, size = 2) + theme_classic()
    gg <- gg + theme(plot.title = element_text(size = 16))
    gg <- gg + theme(axis.text.x = element_text(size = 12))
    gg <- gg + theme(axis.text.y = element_text(size = 12))
    gg <- gg + theme(legend.position = "none")
    gg <- gg + labs(title = plot.title)  + xlab(xlab) + ylab(ylab)
    #print(gg)
    return(gg)
  } #end hazard

}
