library(rstan)
library(parallel)
library(dplyr)

rm(list = ls())
setwd("~/Dropbox/GitHub/glmm")

load("~/Dropbox/GitHub/glmm/data/mdata.rdt")
load("~/Dropbox/GitHub/glmm/data/kinship.rdt")

load("~/Dropbox/GitHub/glmm/Apoe/apoe.rdt")
load("/data/xwang/glmm/prune.rdt")
geno <- apoe[rownames(apoe) %in% prune$V1, ]

source("~/Dropbox/GitHub/glmm/R/waic.R")  # WAIC

chr = "chr19"

cov <- mdata[c("Age", "Sex")]
cov$Sex <- as.integer(cov$Sex) - 1
Ad <- as.numeric(mdata$AD2)
g0 <- rep(0, length(Ad))

Sigma <- kinship[[paste0("no_", chr)]]
Sigma[Sigma < 0] <- 0
L <- t(chol(Sigma))

model <- stan_model("~/Dropbox/GitHub/glmm/Stan/glmm.stan")
dat = list(N = 576, K = 4, D = 2, cov = cov, L = L, g = g0, Ad = Ad)

myGWAS_sampling <- function(x) {
  y1 = list()
  N1 = nrow(x)

  for (i in 1:N1) {
    dat1 = within(dat, {g = x[i, ]})
    fit1 = sampling(model, data = dat1, chain = 2, iter = 400, warmup = 200)
    y1$summary[[i]] = summary(fit1, pars = c("p", "beta", "sigma", "c"))$summary
    y1$waic[[i]] = Waic(fit1)$WAIC
  }

  names(y1$summary) = names(y1$waic) = rownames(x)
  return(y1)
}

n.core <- 20 
u.core <- round(nrow(geno) / n.core)
idx1 <- ((1:n.core) - 1) * u.core + 1 
idx2 <- c((1:(n.core - 1)) * u.core, nrow(geno))
genoList <- mclapply(1:15, function(x) geno[idx1[x]:idx2[x], ], mc.preschedule = F)

y <- mclapply(genoList, myGWAS_sampling, mc.cores = n.core)
save(y, file = "/data/xwang/glmm/sampling_Apoe.rdt")

