library(rstan)
library(parallel)
library(dplyr)

rm(list = ls())

arg <- commandArgs(TRUE)
chr <- paste("chr", arg, sep = "")

load("~/Dropbox/GitHub/glmm/data/mdata.rdt")
load("~/Dropbox/GitHub/glmm/data/kinship.rdt")
source("~/Dropbox/GitHub/glmm/R/waic.R")  # WAIC

load(paste0("/data/xwang/glmm/geno/", chr, ".rdt"))
load(paste0("/data/xwang/glmm/prior/", chr, ".rdt"))

load("/data/xwang/glmm/noK_p0.01.rdt")
noK_v <- rownames(noK_p0.01)[noK_p0.01$P < 0.0001]

igap_p <- pnorm(abs(igap$Beta), sd = igap$SE, lower.tail = F) * 2
igap_v <- rownames(igap)[igap_p < 0.0001]

geno <- geno[rownames(geno) %in% c(igap_v, noK_v), ]

cov <- mdata[c("Age", "Sex")]
cov$Sex <- as.integer(cov$Sex) - 1
Ad <- as.numeric(mdata$AD2)
g0 <- rep(0, length(Ad))

Sigma <- kinship[[paste0("no_", chr)]]
Sigma[Sigma < 0] <- 0
L <- t(chol(Sigma))

model <- stan_model("~/Dropbox/GitHub/glmm/Stan/glmm_prior.stan")
dat0 = list(N = 576, K = 4, D = 2, cov = cov, L = L, g = g0, prior = 0, Ad = Ad)

myGWAS_sampling <- function(geno) {
  y1 = list()
  N1 = nrow(geno)
  vId = rownames(geno)
  p1 = igap[vId, 1] / igap[vId, 2]

  for (i in 1:N1) {
    dat1 = within(dat0, { g = geno[i, ]; prior = p1[i] })
    fit = sampling(model, data = dat1, chain = 3, iter = 400, warmup = 200)
    y1$summary[[i]] = summary(fit, pars = c("p", "beta"))$summary
#   y1$waic[[i]] = Waic(fit)$WAIC
  }

  names(y1$summary) = rownames(geno)
# names(y1$summary) = names(y1$waic) = rownames(x)
  return(y1)
}

n.core <- 60 
u.core <- round(nrow(geno) / n.core)
idx1 <- ((1:n.core) - 1) * u.core + 1 
idx2 <- c((1:(n.core - 1)) * u.core, nrow(geno))
genoList <- mclapply(1:n.core, function(x) geno[idx1[x]:idx2[x], ], mc.preschedule = F)

y <- mclapply(genoList[1:15], myGWAS_sampling, mc.cores = 15)
save(y, file = paste0("/data/xwang/glmm/sampling_prior/", chr, "_p1.rdt"))

y <- mclapply(genoList[16:30], myGWAS_sampling, mc.cores = 15)
save(y, file = paste0("/data/xwang/glmm/sampling_prior/", chr, "_p2.rdt"))

y <- mclapply(genoList[31:45], myGWAS_sampling, mc.cores = 15)
save(y, file = paste0("/data/xwang/glmm/sampling_prior/", chr, "_p3.rdt"))

y <- mclapply(genoList[46:60], myGWAS_sampling, mc.cores = 15)
save(y, file = paste0("/data/xwang/glmm/sampling_prior/", chr, "_p4.rdt"))

