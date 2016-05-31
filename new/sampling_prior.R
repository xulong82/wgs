library(rstan)
library(parallel)
library(dplyr)

rm(list = ls())

chr <- commandArgs(TRUE)
glmm_prior <- stan_model("~/Dropbox/GitHub/wgs/new/glmm_prior.stan")

load("~/Dropbox/GitHub/wgs/new/mdata.rdt")

load("/data/xwang/adsp3/igap/igap.rdt")
igap <- igap_s1[igap_s1$Chromosome == chr, ]
igap$UID <- paste0(igap[, 1], ":", igap[, 2], "_", igap[, 5], "/", igap[, 4])

load(paste0("/data/xwang/adsp3/R/chr", chr, ".rdt"))
geno <- geno[rownames(geno) %in% igap$UID, ]

load("/data/xwang/adsp3/glm.rdt")

v1 <- rownames(glm)[glm$P < 0.0001]
v2 <- igap$UID[igap$Pvalue < 0.0001]

geno <- geno[rownames(geno) %in% c(v1, v2), ]
geno <- geno[, mdata$ADSP.Sample.ID]

load("/data/xwang/adsp3/plink/kin.rdt")
kin <- kin[mdata$ADSP.Sample.ID, mdata$ADSP.Sample.ID]

cov <- mdata[c("Age", "Sex")]
dat <- list(N = 570, K = 4, D = 2, prior = 0, cov = cov)
dat <- within(dat, { Ad = as.numeric(mdata$AD1); L = t(chol(kin)) })

myGWAS_sampling <- function(geno) {
  y1 = list()
  for (i in 1:nrow(geno)) {
    dat1 = within(dat, {g = geno[i, ]})
    fit1 = sampling(glmm_prior, data = dat1, chain = 3, iter = 400, warmup = 200)
    y1[[i]] = summary(fit1, pars = c("p", "beta"))$summary
  }
  names(y1) = rownames(geno)
  return(y1)
}

n.core <- 20 
u.core <- round(nrow(geno) / n.core)
idx1 <- ((1:n.core) - 1) * u.core + 1 
idx2 <- c((1:(n.core - 1)) * u.core, nrow(geno))
genoList <- mclapply(1:n.core, function(x) geno[idx1[x]:idx2[x], ], mc.preschedule = F)

sapply(genoList, nrow)

fit <- mclapply(genoList, myGWAS_sampling, mc.cores = 20)
save(fit, file = paste0("/data/xwang/adsp3/sampling_prior/chr", chr, ".rdt"))

