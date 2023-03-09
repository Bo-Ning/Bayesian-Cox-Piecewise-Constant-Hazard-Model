rm(list = ls())

library(survival)
library(ggplot2)
library(dplyr)
library(ggfortify)

##############
# Bayesian method
lung$time.scale <- lung$time/max(lung$time)
lung$cens <- lung$status - 1
lung$trt <- lung$age

df <- data.frame(time = lung$time.scale, eventtime = lung$time.scale, 
                 event = lung$cens, status = lung$cens, trt = scale(lung$age))

cox <- coxph(Surv(time, status) ~ trt, data = df)
summary(cox)
cox_fit <- survfit(cox)
autoplot(cox_fit)

# run MCMC method
source("BayesCox.R")
source("Lambda.at.time.t.R")
source("CoxReshapeData.R")
source("CoxSamplePosteriorDepGamma.R")
source("CoxSamplePosteriorIndepGamma.R")
source("CoxSurvEval.R")
source("CoxSurvivalFromRegAndCumhaz.R")
source("CumhazEval.R")
source("MCMCDepGammaFirst.R")
source("MCMCDepGammaIntermediate.R")
source("PlotBayesSurv.R")

result <- BayesCox(df = df, prior = "Independent", N = 20000, theta.prop.sd = 0.1)
theta.hat <- result$theta.post.mean
surv.grid <- result$surv.eval.grid
surv.grid.data <- surv.grid * max(lung$time)
# surv.grid.data <- surv.grid
surv.lq <- result$surv.lower.quantile
surv.uq <- result$surv.upper.quantile
surv.pm <- result$surv.post.mean

lines(surv.grid.data, surv.lq, col = "blue")
lines(surv.grid.data, surv.uq, col = "blue")
lines(surv.grid.data, surv.pm, col = "red")
plot(cox_fit$time*max(lung$time), cox_fit$lower, type = "l", ylab = " ", xlab = "Time", main = "")
lines(cox_fit$time*max(lung$time), cox_fit$upper, type = "l")

bayes.data <- data.frame(time = surv.grid.data,
                         bayes.lower = surv.lq, bayes.upper = surv.uq,
                         bayes.mean = surv.pm)
freq.data <- data.frame(surv.grid = cox_fit$time*max(lung$time),
                        freq.lower = cox_fit$lower, freq.upper = cox_fit$upper)
library(ggplot2)
ggplot(bayes.data) + geom_line(aes(y=bayes.mean, x=time), colour = "red", size = 1)+
  geom_ribbon(aes(ymin=bayes.upper, ymax=bayes.lower, x=time, fill = "95% Bayesian credible band"), alpha = 0.5) +
  scale_colour_manual("",values="blue")+
  scale_fill_manual("",values="orange")+ 
  geom_line(aes(y = freq.lower, x=surv.grid, color="95% Kaplan-Meier band"), data = freq.data, colour = "blue", size = 1, 
            linetype = "longdash") +
  geom_line(aes(y = freq.upper, x=surv.grid), data = freq.data, colour = "blue", size = 1, linetype = "longdash") +
  xlab("Time") + ylab(" ") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )
  



###############################################
rm(list = ls())
library(spBayesSurv)
Leuk <- data(LeukSurv)
Leuk <- LeukSurv

Leuk$time.scale <- Leuk$time/max(Leuk$time)
Leuk$trt <- Leuk$age

cox <- coxph(Surv(time, cens) ~ age, data = LeukSurv)
summary(cox)
cox_fit <- survfit(cox)
autoplot(cox_fit)
plot(cox_fit$time, cox_fit$lower, type = "l", xlab = "Time", ylab = " ")
lines(cox_fit$time, cox_fit$upper, type = "l")


df <- data.frame(time = Leuk$time.scale, eventtime = Leuk$time.scale, 
                 event = Leuk$cens, status = Leuk$cens, trt = Leuk$age)

# run MCMC method
source("BayesCox.R")
source("Lambda.at.time.t.R")
source("CoxReshapeData.R")
source("CoxSamplePosteriorDepGamma.R")
source("CoxSamplePosteriorIndepGamma.R")
source("CoxSurvEval.R")
source("CoxSurvivalFromRegAndCumhaz.R")
source("CumhazEval.R")
source("MCMCDepGammaFirst.R")
source("MCMCDepGammaIntermediate.R")
source("PlotBayesSurv.R")

result <- BayesCox(df = df, prior = "Independent", N = 20000, theta.prop.sd = 0.1)
theta.hat <- result$theta.post.mean
surv.grid <- result$surv.eval.grid
surv.grid.data <- surv.grid * max(Leuk$time)
# surv.grid.data <- surv.grid
surv.lq <- result$surv.lower.quantile
surv.uq <- result$surv.upper.quantile
surv.pm <- result$surv.post.mean
lines(surv.grid.data, surv.lq, col = "blue")
lines(surv.grid.data, surv.uq, col = "blue")
lines(surv.grid.data, surv.pm, col = "red")

