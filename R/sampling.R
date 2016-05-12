library(rstan)
library(parallel)
library(dplyr)

rm(list = ls())

arg <- commandArgs(TRUE)
chr <- paste("chr", arg, sep = "")

load("~/Dropbox/GitHub/glmm/data/mdata.rdt")
load("~/Dropbox/GitHub/glmm/data/kinship.rdt")
source("~/Dropbox/GitHub/glmm/R/waic.R")  # WAIC

load(paste0("/data/xwang/Adsp/golden/geno_", chr, ".rdt"))

load("/data/xwang/glmm/noK_p0.01.rdt")
noK_v <- rownames(noK_p0.01)[noK_p0.01$P < 0.0001]

geno <- geno[rownames(geno) %in% noK_v, ]

cov <- mdata[c("Age", "Sex")]
cov$Sex <- as.integer(cov$Sex) - 1
Ad <- as.numeric(mdata$AD2)
g0 <- rep(0, length(Ad))

Sigma <- kinship[[paste0("no_", chr)]]
Sigma[Sigma < 0] <- 0
L <- t(chol(Sigma))

model <- stan_model("~/Dropbox/GitHub/glmm/Stan/glmm.stan")
dat0 = list(N = 576, K = 4, D = 2, cov = cov, L = L, g = g0, Ad = Ad)

myGWAS_sampling <- function(geno) {
  y1 = list()
  N1 = nrow(geno)

  for (i in 1:N1) {
    dat1 = within(dat0, {g = geno[i, ]})
    fit = sampling(model, data = dat1, chain = 3, iter = 400, warmup = 200)
    y1$summary[[i]] = summary(fit, pars = c("p", "beta"))$summary
#   y1$waic[[i]] = Waic(fit)$WAIC
  }

  names(y1$summary) = rownames(geno)
# names(y1$summary) = names(y1$waic) = rownames(x)
  return(y1)
}

n.core <- 40 
u.core <- round(nrow(geno) / n.core)
idx1 <- ((1:n.core) - 1) * u.core + 1 
idx2 <- c((1:(n.core - 1)) * u.core, nrow(geno))
genoList <- mclapply(1:40, function(x) geno[idx1[x]:idx2[x], ], mc.preschedule = F)

y <- mclapply(genoList[1:20], myGWAS_sampling, mc.cores = 20)
save(y, file = paste0("/data/xwang/glmm/sampling/", chr, "_p1.rdt"))

y <- mclapply(genoList[21:40], myGWAS_sampling, mc.cores = 20)
save(y, file = paste0("/data/xwang/glmm/sampling/", chr, "_p2.rdt"))

