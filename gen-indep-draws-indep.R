# run 1000 MCMC draws
setwd("~/Cox-simulation")
num.obs <- 2000
sim.num <- 1000 # simulation 1000 datasets
cores <- 20

# load libraries
library(parallel)
library(simsurv)
library(ggplot2)
library(BayesSurvival)
library(survival)

# load functions
source("BayesCox.R")
source("CoxReshapeData.R")
source("CoxSamplePosteriorDepGamma.R")
source("CoxSamplePosteriorIndepGamma.R")
source("CoxSurvEval.R")
source("CoxSurvivalFromRegAndCumhaz.R")
source("CumhazEval.R")
source("MCMCDepGammaFirst.R")
source("MCMCDepGammaIntermediate.R")
source("PlotBayesSurv.R")
source("Lambda.at.time.t.R")

set.seed(20211003)

cond.hazard.true <- function(t, x, betas, ...){(6*((t + 0.05)^3 - 2*(t + 0.05)^2 + t + 0.05) + 0.7)*exp(x*betas)}
covs <- data.frame(id = 1:num.obs, trt = stats::rnorm(num.obs, 0, 1))
df.tmp <- simsurv(hazard = cond.hazard.true, betas = c(trt = -0.5), x = covs, maxt = 1)
cens <- runif(num.obs, min = 0, max = 1)
df <- data.frame(eventtime = pmin(df.tmp$eventtime, cens),
                 status = as.numeric(df.tmp$eventtime <= cens))
df$time <- df$eventtime
df$event <- df$status
df$trt <- covs$trt

simulation <- function(num.obs, prior, N, df) {
  # prepare for MCMC
  p <- 1
  model <- coxph(Surv(eventtime, status) ~ trt, data=df)
  theta.int <- model$coefficients # frequentist estimaton
  exp.reg <- c(exp(df$trt * theta.int))
  bsline.cumhaz.coxph <- basehaz(model, centered = F) # true cumulative hazard
  K = ceiling((dim(df)[1] / log(dim(df)[1]))^(1/2))
  interval <- c(0,cumsum(rep(1/K, K)))
  Lambda.int <- sapply(1:K, function(k) {
    if (k == 1) {
      index <- bsline.cumhaz.coxph$time >= interval[k] & bsline.cumhaz.coxph$time <= interval[k+1]
      hzd.val <- bsline.cumhaz.coxph$hazard[index]
      time.val <- bsline.cumhaz.coxph$time[index]
      length.diff <- time.val - c(0, time.val[1:(sum(index)-1)])
      val <- sum(hzd.val * length.diff)/sum(length.diff)
    } else {
      index <- bsline.cumhaz.coxph$time >= interval[k] & bsline.cumhaz.coxph$time <= interval[k+1]
      hzd.val <- bsline.cumhaz.coxph$hazard[index]
      time.val <- bsline.cumhaz.coxph$time[index]
      length.diff <- time.val - c(bsline.cumhaz.coxph$time[which(index == T)[1] - 1], time.val[1:(sum(index)-1)])
      val <- sum(hzd.val * length.diff)/sum(length.diff)
    }
  })
  lambda.int <- Lambda.int - c(0,Lambda.int[1:(K-1)])
  time.max = max(df$time)
  data.summary <- CoxReshapeData(df, time, event, exp.reg, lambda.int, K, time.max)
  fail <- data.summary$failures
  exp <- data.summary$exposures
  # run MCMC
  post.samples <- SamplePosteriorIndepGamma(fail, exp, df, theta.int, p,
                                            N = N, 
                                            alpha.indep = 1, beta.indep = 1, 
                                            theta.prop.sd = 0.2)
  
  # collect result
  bsline.haz.samples <- post.samples$lambda.samples[N, ]
  theta.samples <- post.samples$theta.samples[N]
  return(list(bsline.haz.samples = bsline.haz.samples,
              theta.samples = theta.samples))
}

Results <- mclapply(1:sim.num, simulation, prior = "Inependent", N = 10000,
                    df = df, mc.cores = cores)

theta.draws <- haz.draws <- c()

for (i in 1:sim.num) {
  theta.draws <- cbind(theta.draws, Results[[i]]$theta.samples)
  haz.draws <- cbind(haz.draws, Results[[i]]$bsline.haz.samples)
}

write.csv(theta.draws, file = "theta-draws-ind.csv")
write.csv(haz.draws, file = "haz-draws-ind.csv")

