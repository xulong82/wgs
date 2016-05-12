#  Bayesian inference on ADSP project
#  -------------------------------------------------------------------
#  Random effect: Genotype-based Relatedness

library(rstan)
library(parallel)
library(mail)

rm(list = ls())

# --- user configurations ---
nchains = 2  # number of chains
nwarmup = 3e2  # warm up (discarded)
niter = 6e2  # number of iterations (including warmup)

# --- META data ---
load("~/Dropbox/ADSP/Stan/pdata.rdt")

# --- GWAS
source("~/Dropbox/ADSP/Stan/waic.R")

load("/data/xwang/ADSP/genotype1/autosome.rdt")
geno <- geno[snpInf$CHROM == "chr2", ]

prune.in <- read.delim("/data/xwang/ADSP/kinship/plink.prune.in", 
  stringsAsFactors = F, header = F)$V1
geno <- geno[rownames(geno) %in% prune.in, ]

load("~/Dropbox/ADSP/Stan/Umode.rdt")
U0 <- Umode[[2]]

model_glm <- stan_model(file = "~/Dropbox/ADSP/Stan/glm.stan", model_name='GLMM')

cat(date(), "--- split data for node and core \n")
arg <- commandArgs(TRUE)
fname <- paste("glmPUfit", arg, sep = "")
node <- as.numeric(arg)

n.core <- n.node <- 20
u.node <- round(nrow(geno) / n.node)
idx1 <- ((1:n.node) - 1) * u.node + 1
idx2 <- c((1:(n.node - 1)) * u.node, nrow(geno))
geno <- geno[idx1[node]:idx2[node], ]
u.core <- round(nrow(geno) / n.core)
idx1 <- ((1:n.core) - 1) * u.core + 1
idx2 <- c((1:(n.core - 1)) * u.core, nrow(geno))

genoList <- list()
for (i in 1:n.core)
  genoList[[i]] <- geno[idx1[i]:idx2[i], ]

myGWAS <- function(geno) {
  fit = list()
  for (i in 1:nrow(geno)) {
#   cat(i, "in", nrow(geno), rownames(geno)[i], "\n")
    Snp = geno[i, ]
    fit0 = sampling(model_glm, par = c("pSnp", "logLik"), 
                    data = list(N = N, K = K, Ad = Ad, Age = Age, Sex = Sex, Apoe = Apoe, Snp = Snp, U = U0),
#                   par = c("pAge", "pSex", "pSnp", "pApoe22", "pApoe23", "pApoe24", "pApoe34", "pU", "logLik"), 
                    warmup = nwarmup, iter = niter, chains = nchains, refresh = 1e2)
    fit$waic[[i]] = Waic(fit0)$WAIC
    fit$pSnp[[i]] = extract(fit0, pars = "pSnp")  # mode only?
  }
  names(fit$pSnp) <- names(fit$waic) <- rownames(geno)
  return(fit)
}

startTime = Sys.time()
#---------------------------------------------------------------------
x <- mclapply(genoList, myGWAS, mc.cores = n.core)
assign(fname, x)
filepath <- "/data/xwang/ADSP/Stan/chr2"
save(list = fname, file = paste(paste(filepath, fname, sep = "/"), "rdt", sep = "."))
endTime = Sys.time()

timeTook <- endTime - startTime
if(node == 1) sendmail("xulong.wang@jax.org", "R notice", as.character.Date(timeTook))

#---------------------------------------------------------------------
filepath <- "/data/xwang/ADSP/Stan/chr2"
waic <- NULL
for (i in 1:20) {
  fname = paste("glmPUfit", i, sep = "")
  load(paste(paste(filepath, fname, sep = "/"), "rdt", sep = "."))
  x = get(fname)
  for (j in 1:20) {
    waic = c(waic, x[[i]]$waic)
  }
}
load("~/Dropbox/ADSP/GWAS/snpInf.rdt")
pos <- snpInf[match(names(waic), snpInf$ID), "ID"]
save(waic, file = "~/Dropbox/ADSP/Stan/waic2.rdt")
