setwd("~/chzf/")
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
library(mvtnorm)

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
simulation <- function(num.obs, prior, N, p) {
  cond.hazard.true <- function(t, x, betas, ...) {(sin((t+0.05)*2*pi)*0.8 + (t + 0.05)^4 - 1.8*(t + 0.05)^2 + 2) *     
    exp(as.matrix(x) %*% betas)}
  covs <- data.frame(trt = matrix(stats::rnorm(num.obs*p, 0, 1), num.obs, p))
  df <- simsurv(hazard = cond.hazard.true, 
                betas = c(trt = c(-0.3, -0.1, 0, 0.1, 0.3)), 
                x = covs, maxt = 1)
  df$time <- df$eventtime
  df$event <- df$status
  df$trt <- as.matrix(covs)
  result <- BayesCox(df = df, prior = prior, N = N)
  theta.qtl <- result$theta.qtl
  bs.cumhaz.qtl <- result$bs.cumhaz.qtl
  cond.cumhaz.qtl <- result$cond.cumhaz.qtl
  bs.surv.qtl <- result$bs.surv.qtl
  cond.surv.qtl <- result$cond.surv.qtl
  if (p > 1) {
    z.bar <- colMeans(df$trt)
  } else {
    z.bar <- mean(df$trt)
  }
  return(list(theta.qtl = theta.qtl, 
              bs.cumhaz.qtl = bs.cumhaz.qtl,
              cond.cumhaz.qtl = cond.cumhaz.qtl,
              bs.surv.qtl = bs.surv.qtl,
              cond.surv.qtl = cond.surv.qtl,
              time.max = result$time.max, 
              z.bar = z.bar)
  )
}

# parallel computing
p <- 5
Results <- mclapply(rep(num.obs, sim.num), simulation, prior = "Dependent", 
                    N = 10000, p = p, mc.cores = cores)

theta.qlt.sim <- bs.cumhaz.qtl.sim <- cond.cumhaz.qtl.sim <- 
  bs.surv.qtl.sim <- cond.surv.qtl.sim <- time.max.sim <- z.bar.sim <- c()

for (i in (1:sim.num)) {
  theta.qlt.sim <- rbind(theta.qlt.sim, Results[[i]]$theta.qtl)
  bs.cumhaz.qtl.sim <- rbind(bs.cumhaz.qtl.sim, Results[[i]]$bs.cumhaz.qtl)
  cond.cumhaz.qtl.sim <- rbind(cond.cumhaz.qtl.sim, Results[[i]]$cond.cumhaz.qtl)
  bs.surv.qtl.sim <- rbind(bs.surv.qtl.sim, Results[[i]]$bs.surv.qtl)
  cond.surv.qtl.sim <- rbind(cond.surv.qtl.sim, Results[[i]]$cond.surv.qtl)
  time.max.sim <- cbind(time.max.sim, Results[[i]]$time.max)
  z.bar.sim <- cbind(z.bar.sim, Results[[i]]$z.bar)
}

# theta.est <- matrix(theta.qlt.sim$theta.pm, p, sim.num)
# theta.est.lq <- matrix(theta.qlt.sim$theta.lq, p, sim.num)
# theta.est.uq <- matrix(theta.qlt.sim$theta.uq, p, sim.num)
# rowMeans(theta.est)
# rowMeans(theta.est.lq)
# rowMeans(theta.est.uq)

write.csv(theta.qlt.sim, file = "depAdm200-theta.qlt.sim.csv")
write.csv(bs.cumhaz.qtl.sim, file = "depAdm200-bs.cumhaz.qtl.sim.csv")
write.csv(cond.cumhaz.qtl.sim, file = "depAdm200-cond.cumhaz.qtl.sim.csv")
write.csv(bs.surv.qtl.sim, file = "depAdm200-bs.surv.qtl.sim.csv")
write.csv(cond.surv.qtl.sim, file = "depAdm200-cond.surv.qtl.sim.csv")
write.csv(time.max.sim, file = "depAdm200-time.max.sim.csv")
write.csv(z.bar.sim, file = "depAdm200-z.bar.sim.csv")


