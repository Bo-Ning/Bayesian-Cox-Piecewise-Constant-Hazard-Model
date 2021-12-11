# read file
rm(list = ls())
url <- "depAdm200"
num.hist <- 12

setwd(paste("~/Desktop/research/BayesianCoxModel/Simulation/simulation-result/result-w-condSurv-beta-quarter/", url, sep = ""))
# setwd("~/Desktop/Cox-simulation/")
surv.qlt <- read.csv(paste(url,"-bs.surv.qtl.sim.csv",sep=""), header = T)
bs.surv.pm <- matrix(surv.qlt$bs.surv.pm, num.hist*10, 1000)
bs.surv.lb <- matrix(surv.qlt$bs.surv.lb, num.hist*10, 1000)
bs.surv.ub <- matrix(surv.qlt$bs.surv.ub, num.hist*10, 1000)
surv.grid <- matrix(surv.qlt$surv.eval.grid, num.hist*10, 1000)
# true survival function 
hazard.true <- function(t, x, betas, ...){(6*((t + 0.05)^3 - 2*(t + 0.05)^2 + t + 0.05) + 0.7)}
cumhaz.true <- Vectorize( function(t){integrate(hazard.true, 0, t)$value} )
surv.true <- function(t){exp(-integrate(hazard.true, 0, t)$value)}


# calculate confidence band for Bayesian method
coverage <- 0
area <- c()
for (i in 1:1000) {
  surv.true.val <- sapply(surv.grid[, i], surv.true)
  lower <- pmax(0, bs.surv.lb[, i])
  upper <- pmin(1, bs.surv.ub[, i])
  diff.lq <- surv.true.val - lower
  diff.uq <- upper - surv.true.val 
  area <- c(area, sum(upper - lower)/(num.hist*10))
  if ((min(diff.lq)) >= 0 & (min(diff.uq)) >= 0) {
    coverage <- coverage+1
    # cat("i= ", i, "coverage = ", coverage, "\n")
  }
}

cat("coverage = ", coverage, ", area = ", mean(area))


cond.surv.qlt <- read.csv(paste(url,"-cond.surv.qtl.sim.csv",sep=""), header = T)
cond.surv.pm <- matrix(surv.qlt$bs.surv.pm, num.hist*10, 1000)
cond.surv.lb <- matrix(surv.qlt$bs.surv.lb, num.hist*10, 1000)
cond.surv.ub <- matrix(surv.qlt$bs.surv.ub, num.hist*10, 1000)
surv.grid <- matrix(surv.qlt$surv.eval.grid, num.hist*10, 1000)
z <- read.csv(paste(url,"-z.bar.sim.csv", sep = ""), header = T)
# calculate confidence band for Bayesian method
coverage.z <- 0
area.z <- c()
for (i in 1:1000) {
  cumhaz.true.val <- sapply(surv.grid[, i], cumhaz.true)
  cond.surv.true.val <- exp(-cumhaz.true.val*exp(-0.5*z[[i+1]]))
  lower <- pmax(0, cond.surv.lb[,i])
  upper <- pmin(1, cond.surv.ub[,i])
  diff.lq <- cond.surv.true.val - lower
  diff.uq <- upper - cond.surv.true.val 
  area.z <- c(area.z, sum(upper - lower)/(num.hist*10))
  if ((min(diff.lq)) >= 0 & (min(diff.uq)) >= 0) {
    coverage.z <- coverage.z + 1
    # cat("i= ", i, "coverage = ", coverage, "\n")
  }
}

cat("coverage.z = ", coverage.z, ", area.z = ", mean(area.z))


plot(surv.grid[, i], lower, type = "l", ylim = c(0, 1))
lines(surv.grid[, i], upper, type = "l")
lines(surv.grid[, i], surv.true.val, type = "l", col = "red")

