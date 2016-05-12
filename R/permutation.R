library(rstan)
library(parallel)
library(dplyr)

rm(list = ls())
setwd("~/Dropbox/GitHub/Adsp")

load("./data/sampling.rdt")

arg <- commandArgs(TRUE) 
chr <- paste("chr", arg, sep = "")

# glm_stan <- stan_model("./GLM/glm(noK).stan")
glm_stan <- stan_model("./GLM/glm(151024).stan")

cov <- data$mdata[c("Age", "Sex", "Apoe2", "Apoe4")]
cov$Var <- rep(0, nrow(cov))
cov$Sex <- as.integer(cov$Sex) - 1
Ad <- as.numeric(data$mdata$AD2)

init = data$init

Sigma <- data$kinship[[paste0("no_chr", arg)]]
Sigma[Sigma < 0] <- 0
L <- t(chol(Sigma))

dat = list(N = 576, K = 4, D = 5, cov = cov, Ad = Ad)

myGWAS_optimizing_permutation <- function(idx) {

  N1 = 1e5
  maxLogLik = -1e3
  geno0 = geno[sample(1:nrow(geno), N1), ]

  for (i in 1:N1) {
    dat$cov$Var = sample(geno0[i, ])
    fit = optimizing(glm_stan, algorithm = "LBFGS", init = init, data = dat)
    loglik = sum(fit$par[paste0("lp[", 1:576, "]")])

    if (loglik > maxLogLik) maxLogLik = loglik
  }
  return(maxLogLik)

}

setwd("/data/xwang/Adsp")
load(paste0("./golden/geno_", chr, ".rdt"))

y <- mclapply(1:3e1, myGWAS_optimizing_permutation, mc.cores = 15)
assign(chr, unlist(y))

# save(list = chr, file = paste0("./permutation/R1/", chr, ".rdt"))
# save(list = chr, file = paste0("./permutation/R2/", chr, ".rdt"))
# save(list = chr, file = paste0("./permutation/R3/", chr, ".rdt"))
  save(list = chr, file = paste0("./permutation/R4/", chr, ".rdt"))

###

# x = replicate(100, {
#   par0 = optimizing(glm_stan, verbose = FALSE, algorithm = "LBFGS", data = dat)$par
#   sum(par0[paste("lp[", 1:576, "]", sep = "")])
# })
# 
# init <- list(c =  par0[c("c[1]", "c[2]", "c[3]")],
#              beta = par0[c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]")])

