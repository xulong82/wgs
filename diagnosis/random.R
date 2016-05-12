library(rstan)
library(parallel)
library(mail)

rm(list = ls())

load("~/Dropbox/ADSP/Stan/pdata.rdt")  # phenotype
load("~/Dropbox/ADSP/Stan/geno_rs12721046.rdt")  # APOE locus
geno <- geno[131:530, ]  # 400 SNP with KS random effect inference
x = snpInf[match(rownames(geno), snpInf$ID), "POS"] * 1e-3

load("~/Dropbox/ADSP/kinship/kinship.rdt")  # kinship matrix
KS <- kinship$chr[[19]]
KS[KS < 0] <- 0

# load("~/Dropbox/ADSP/kinship/KS_doqtl.rdt")
# KS <- kin2

# model_glm <- stan_model(file = "~/Dropbox/ADSP/Stan/random/glmFixApoeNull.stan", model_name='GLMM')
# model_glm <- stan_model(file = "~/Dropbox/ADSP/Stan/random/glmKS1ApoeNull.stan", model_name='GLMM')
# model_glm <- stan_model(file = "~/Dropbox/ADSP/Stan/random/glmKS2ApoeNull.stan", model_name='GLMM')
  model_glm <- stan_model(file = "~/Dropbox/ADSP/Stan/random/glmKS3ApoeNull.stan", model_name='GLMM')
# model_glm <- stan_model(file = "~/Dropbox/ADSP/Stan/random/glmPUApoeNull.stan", model_name='GLMM')

# --- Maximum Likelihood Optimizing ---
fit0 <- optimizing(model_glm, verbose = FALSE, algorithm = "LBFGS",
                   data = list(N = N, K = K, Ad = Ad, Age = Age, Sex = Sex, Snp = rep(0, N), KS = KS))
genInitList <- function() {
  list(z = fit0$par[paste(paste("z[", 1:576, sep = ""), "]", sep = "")], 
       sigmaU = fit0$par["sigmaU"], 
       alpha = fit0$par["alpha"], 
       pAge = fit0$par["pAge"], 
       pSex = fit0$par["pSex"], 
       cut1 =  fit0$par[c("cut1[1]", "cut1[2]", "cut1[3]")],
       pSnp = 0) 
}
save(fit0, genInitList, file = "~/Dropbox/ADSP/Stan/random/genInitList.rdt")

loglik <- NULL
mypars <- paste(paste("logLik[", 1:576, sep = ""), "]", sep = "")
load("~/Dropbox/ADSP/Stan/random/genInitList.rdt")

for (i in 1:nrow(geno)) {
  if (i %% 10 == 0) cat(i, "in", nrow(geno), "\n")
  geno1 <- geno[i, ]
  fit1 <- optimizing(model_glm, verbose = FALSE, algorithm = "LBFGS",  init = genInitList,
                     data = list(N = N, K = K, Ad = Ad, Age = Age, Sex = Sex, Snp = geno1, KS = KS))
  
  loglik1 = 0
  for (j in 1:576) loglik1 = loglik1 + fit1$par[[mypars[j]]] 
  loglik <- c(loglik, loglik1)
}

plot(x, 2* (loglik - min(loglik)))

pdf("~/Dropbox/ADSP/Stan/random/glmKS3Apoe.pdf", height = 4)
plot(x, 2* (loglik - min(loglik)),
     main = "APOE", xlab = "", ylab = "LOD",
     ylim = c(0, 27),
     type = "p", pch = 19, cex = 1, col = "lightblue")
dev.off()

sigmaU <- NULL
for (i in 1:length(fit)) {
  sigmaU <- c(sigmaU, fit[[i]]$par["sigmaU"])
}

# --- MCMC Sampling ---
arg <- commandArgs(TRUE)

node <- as.numeric(arg)

n.core <- n.node <- 20
geno <- geno[((node - 1) * 20 + 1):(node * 20), ]

genoList <- list()
for (i in 1:n.core) genoList[[i]] <- geno[i, ]

source("~/Dropbox/ADSP/Stan/waic.R")
load("~/Dropbox/ADSP/Stan/random/genInitList.rdt")

myGWAS <- function(geno) {
  fit0 = sampling(model_glm, init = genInitList,  # par = c("pSnp", "pU", "logLik"), 
#                 data = list(N = N, K = K, Ad = Ad, Age = Age, Sex = Sex, Snp = geno),
                  data = list(N = N, K = K, Ad = Ad, Age = Age, Sex = Sex, Snp = geno, KS = KS),
#                 data = list(N = N, K = K, Ad = Ad, Age = Age, Sex = Sex, Snp = geno, U = Umode),
                  warmup = 3e2, iter = 6e2, chains = 2, refresh = 1e2)
  return(fit0)
}

startTime = Sys.time()
# ---------------------------------------------------------------------
x <- mclapply(genoList, myGWAS, mc.cores = n.core)
filepath <- "/data/xwang/ADSP/Stan/random"
fname <- paste("glmKS3ApoeNull_fit", arg, sep = "")  # random's full estimate with kinship
assign(fname, x)
# fname <- paste("glmPUApoeNull_fit", arg, sep = "")  # random's point estimate

save(list = fname, file = paste(paste(filepath, fname, sep = "/"), "rdt", sep = "."))
endTime = Sys.time()

timeTook <- endTime - startTime
if(node == n.node) sendmail("xulong.wang@jax.org", "R notice", as.character.Date(timeTook))

# --- extract posterior sampling
par1 <- c("pAge", "pSex", "pApoe22", "pApoe23", "pApoe24", "pApoe34")
par2 <- c(par1, "cut[1]", "cut[2]", "cut[3]")
par3 <- paste(paste("u[", 1:576, sep = ""), "]", sep = "")
  
pSnp <- rep(NULL, 6e2)
waic <- NULL

filepath <- "/data/xwang/ADSP/Stan/random"
for (i in 1:20) {
  fname = paste("glmKS3ApoeNull_fit", i, sep = "")
# fname = paste("glmPUApoeNull_fit", i, sep = "")
  cat(fname, "\n")
  load(paste(paste(filepath, fname, sep = "/"), "rdt", sep = "."))
  fit0 <- get(fname)
  for (j in 1:20) {
    pSnp <- cbind(pSnp, unlist(extract(fit0[[j]], pars = "pSnp")))
    waic <- c(waic, Waic(fit0[[j]])$WAIC)
  }
}

save(pSnp, waic, file = "~/Dropbox/ADSP/Stan/random/glmKS3ApoeNull.rdt")
save(pSnp, waic, file = "~/Dropbox/ADSP/Stan/random/glmPUApoeNull.rdt")

load(file = "~/Dropbox/ADSP/Stan/random/glmKS3ApoeNull.rdt")
plot(x, waic)

pdf("~/Dropbox/ADSP/Stan/random/glmKS3ApoeNull_sampling.pdf", height = 4)
plot(x, -waic + max(waic),
     main = "APOE", xlab = "", ylab = "WAIC",
#    ylim = c(0, 35),
     type = "p", pch = 19, cex = 1, col = "lightblue")
dev.off()
