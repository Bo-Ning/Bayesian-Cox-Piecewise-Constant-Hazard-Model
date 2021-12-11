setwd("~/Cox-simulation/")
rm(list = ls())

num.obs <- 200
sim.num <- 1000 # simulation 1000 datasets
cores <- 10

# load libraries
library(parallel)
library(simsurv)
library(ggplot2)
library(BayesSurvival)
library(survival)

# load functions
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

set.seed(20210010)
simulation <- function(num.obs, prior, N) {
  cond.hazard.true <- function(t, x, betas, ...){(6*((t + 0.05)^3 - 2*(t + 0.05)^2 + t + 0.05) + 0.7)*exp(x*betas)}
  covs <- data.frame(id = 1:num.obs, trt = stats::rnorm(num.obs, 0, 1))
  df <- simsurv(hazard = cond.hazard.true, betas = c(trt = -0.5), x = covs, maxt = 1)
  df$time <- df$eventtime
  df$event <- df$status
  df$trt <- covs$trt
  result <- BayesCox(df = df, prior = prior, N = N)
  return(list(theta.post.mean = result$theta.post.mean,
              theta.radius = result$theta.radius,
              theta.lower.quantile = result$theta.lower.quantile, 
              theta.upper.quantile = result$theta.upper.quantile,
              baseline.haz.post.mean = result$baseline.haz.post.mean,
              baseline.cumhaz.post.mean = result$baseline.cumhaz.post.mean,
              baseline.cumhaz.radius = result$baseline.cumhaz.radius,
              baseline.cumhaz.lower.quantile = result$baseline.cumhaz.lower.quantile, 
              baseline.cumhaz.upper.quantile = result$baseline.cumhaz.upper.quantile, 
              surv.post.mean = result$surv.post.mean,
              surv.radius = result$surv.radius,
              surv.lower.quantile = result$surv.lower.quantile,
              surv.upper.quantile = result$surv.upper.quantile,
              surv.eval.grid = result$surv.eval.grid,
              time.max = result$time.max)
  )
}

# parallel computing
Results <- mclapply(rep(num.obs, sim.num), simulation, prior = "Dependent", N = 10000, mc.cores = cores)

theta.post.mean.sim <- theta.radius.sim <- theta.lower.quantile.sim <- 
theta.upper.quantile.sim <- baseline.haz.post.mean.sim <- 
baseline.cumhaz.post.mean.sim <- baseline.cumhaz.radius.sim <-
baseline.cumhaz.lower.quantile.sim <- baseline.cumhaz.upper.quantile.sim <- 
surv.post.mean.sim <- surv.radius.sim <- surv.lower.quantile.sim <-
surv.upper.quantile.sim <- surv.eval.grid.sim <- time.max.sim <- c()

for (i in 1:sim.num) {
  theta.post.mean.sim <- cbind(theta.post.mean.sim, Results[[i]]$theta.post.mean)
  theta.radius.sim <- cbind(theta.radius.sim, Results[[i]]$theta.radius)
  theta.lower.quantile.sim <- cbind(theta.lower.quantile.sim, Results[[i]]$theta.lower.quantile)
  theta.upper.quantile.sim <- cbind(theta.upper.quantile.sim, Results[[i]]$theta.upper.quantile)
  baseline.haz.post.mean.sim <- cbind(baseline.haz.post.mean.sim, Results[[i]]$baseline.haz.post.mean)
  baseline.cumhaz.post.mean.sim <- cbind(baseline.cumhaz.post.mean.sim, Results[[i]]$baseline.cumhaz.post.mean)
  baseline.cumhaz.radius.sim <- cbind(baseline.cumhaz.radius.sim, Results[[i]]$baseline.cumhaz.radius)
  baseline.cumhaz.lower.quantile.sim <- cbind(baseline.cumhaz.lower.quantile.sim, Results[[i]]$baseline.cumhaz.lower.quantile)
  baseline.cumhaz.upper.quantile.sim <- cbind(baseline.cumhaz.upper.quantile.sim, Results[[i]]$baseline.cumhaz.upper.quantile)
  surv.post.mean.sim <- cbind(surv.post.mean.sim, Results[[i]]$surv.post.mean)
  surv.radius.sim <- cbind(surv.radius.sim, Results[[i]]$surv.radius)
  surv.lower.quantile.sim <- cbind(surv.lower.quantile.sim, Results[[i]]$surv.lower.quantile)
  surv.upper.quantile.sim <- cbind(surv.upper.quantile.sim, Results[[i]]$surv.upper.quantile)
  surv.eval.grid.sim <- cbind(surv.eval.grid.sim, Results[[i]]$surv.eval.grid)
  time.max.sim <- cbind(time.max.sim, Results[[i]]$time.max)
}

write.csv(theta.post.mean.sim, file = "depAdm200-theta.post.mean.sim.csv")
write.csv(theta.radius.sim, file = "depAdm200-theta.radius.sim.csv")
write.csv(theta.lower.quantile.sim, file = "depAdm200-theta.lower.quantile.sim.csv")
write.csv(theta.upper.quantile.sim, file = "depAdm200-theta.upper.quantile.sim.csv")
write.csv(baseline.haz.post.mean.sim, file = "depAdm200-baseline.haz.post.mean.sim.csv")
write.csv(baseline.cumhaz.post.mean.sim, file = "depAdm200-baseline.cumhaz.post.mean.sim.csv")
write.csv(baseline.cumhaz.radius.sim, file = "depAdm200-baseline.cumhaz.radius.sim.csv")
write.csv(baseline.cumhaz.lower.quantile.sim, file = "depAdm200-baseline.cumhaz.lower.quantile.sim.csv")
write.csv(baseline.cumhaz.upper.quantile.sim, file = "depAdm200-baseline.cumhaz.upper.quantile.sim.csv")
write.csv(surv.post.mean.sim, file = "depAdm200-surv.post.mean.sim.csv")
write.csv(surv.radius.sim, file = "depAdm200-surv.radius.sim.csv")
write.csv(surv.lower.quantile.sim, file = "depAdm200-surv.lower.quantile.sim.csv")
write.csv(surv.upper.quantile.sim, file = "depAdm200-surv.upper.quantile.sim.csv")
write.csv(surv.eval.grid.sim, file = "depAdm200-surv.eval.grid.sim.csv")
write.csv(time.max.sim, file = "depAdm200-time.max.sim.csv")
