#  Bayesian inference on ADSP project
#  -------------------------------------------------------------------
#  Random effect: Genotype-based Relatedness

library(rstan)
library(parallel)
library(mail)

rm(list = ls())

#--- user configurations ---
sendEmail = T
nchains =  3  # number of chains
nwarmup = 1e3  # warm up (discarded)
niter =  2e3  # number of iterations (including warmup)

#--- META data ---
load("~/Dropbox/ADSP/R/mdata.rdt")
data1 <- mdata[, c("AD1", "Sex", "Age")] # Sex: 0-M, 1-F
data1$AD2 <- factor(data1$AD1, levels = c("0", "0.25", "0.5", "1"), ordered = T) 

N <- nrow(data1)  # sample number
K <- length(unique(data1$AD1))  # AD categories
Ad <- as.numeric(data1$AD2)
Age <- data1$Age
Sex <- as.numeric(data1$Sex)

load("~/Dropbox/ADSP/R/kinship2.rdt")
KS <- kin2

#--- user defined functions ---
Waic <- function(stanfit) {  #--- WAIC computation
  loglik <- extract (stanfit, "logLik")$logLik
  n1 <- nrow(loglik)  # number of samples
  n2 <- ncol(loglik)  # number of data points
  vars <- colSums((loglik - matrix(colMeans(loglik), n1, n2, byrow = T))^2) / (n1 - 1)
  pwaic <- sum(vars)  # effective parameter number
  lpd <- sum(log(colMeans(exp(loglik))))  # log pointwise predictive density
  waic <- -2 * (lpd - pwaic)
  return(list(WAIC = waic, Pwaic = pwaic, LPD = lpd))
}

genInitList <- function() {
  list(z = rep(0, N), sigmaU = 1, pAge = 0, pSex = 0, pSnp = 0, cut = c(1:3))
}

dataList1 <- c("N", "K", "L", "Ad", "Age", "Sex", "Family")

# --- Null ---
# model_glmNull <- stan_model(file = "~/Dropbox/ADSP/bayes/glmNull.stan",
#                             model_name='GLMM (Null)', verbose = FALSE)
# fit <- sampling(model_glmNull, init = genInitList, data = dataList1,
#                 par = c("pAge", "pSex", "cut", "sigmaU"), 
#                 warmup = nwarmup, iter = niter, chains = nchains, refresh = 1e1)
# waic_null <- Waic(fit)$WAIC
# save(fit, waic_null, file = "/data/xwang/ADSP/STAN/glmNull.rdt")

# -- GWAS ---
model_glm <- stan_model(file = "~/Dropbox/ADSP/Stan/glm.stan", 
                        model_name='GLMM', verbose = FALSE)

load("~/Dropbox/ADSP/Stan/genoSNPrs34572242.rdt")
# geno <- geno[1:400, ]  # test

cat(date(), "--- split data for node and core \n")
arg <- commandArgs(TRUE)
fname <- paste("waic", arg, sep = "")
node <- as.numeric(arg)

n.node <- 30

unit.node <- round(nrow(geno) / n.node)
start <- ((1:n.node) - 1) * unit.node + 1
end <- (1:n.node) * unit.node
end[n.node] <- nrow(geno)
geno <- geno[start[node]:end[node], ]

#---------------------------------------------------------------------
waic <- rep(0, nrow(geno))
names(waic) <- rownames(geno)

startTime = Sys.time()
cat(date(), "--- start \n")

for (i in 1:nrow(geno)) {
  cat(i, "in", nrow(geno), rownames(geno)[i], "\n")
  Snp = geno[i, ]
  dataList2 <- c(dataList1, "Snp")
  
  fit1 <- sampling(model_glm, data = dataList2, init = genInitList, iter = 1, chains = 1)
  fit2 <- mclapply(1:nchains, mc.cores = nchains, FUN = function (chain) {
	    stan(fit = fit1, data = dataList2, warmup = nwarmup, iter = niter, 
            par = c("pAge", "pSex", "pSnp", "cut", "sigmaU", "logLik"), 
	    chains = 1, seed = i, chain_id = chain, refresh = 0)})
  fit <- sflist2stanfit(fit2)

# fit <- sampling(model_glm, init = genInitList, 
#                 data = list(N = N, K = K, KS = KS, Ad = Ad, Age = Age, Sex = Sex, Snp = Snp),
#                 par = c("pAge", "pSex", "pSnp", "cut", "sigmaU", "logLik"), 
#                 warmup = nwarmup, iter = niter, chains = nchains, refresh = 1e2)

  waic[i] <- Waic(fit)$WAIC
}

#---------------------------------------------------------------------
assign(fname, waic)
filepath <- "~/Dropbox/ADSP/Stan"
save(list = fname, file = paste(paste(filepath, fname, sep = "/"), "rdt", sep = "."))

endTime = Sys.time()
timeTook = endTime - startTime

# if (sendEmail) sendmail("xulong.wang@jax.org", "R notice", as.character.Date(timeTook))

