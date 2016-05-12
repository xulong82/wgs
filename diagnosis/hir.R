#  Bayesian inference on ADSP project
#  -------------------------------------------------------------------
#  Random effect: Family-based Relatedness

library(rstan)
library(parallel)

rm(list = ls())

#--- META data
load("~/Dropbox/ADSP/R/mdata.rdt")
data1 <- mdata[, c("AD1", "Sex", "Age", "APOE", "Family.ID")] # Sex: 0-M, 1-F
data1$AD2 <- factor(data1$AD1, levels = c("0", "0.25", "0.5", "1"), ordered = T) 

N <- nrow(data1)  # sample number
K <- length(unique(data1$AD1))  # AD categories
L <- length(unique(data1$Family.ID))  # family number
Ad <- as.numeric(data1$AD2)
Age = data1$Age
Sex = as.numeric(data1$Sex)
Family <- as.numeric(data1$Family.ID)  # 1 -> 111

#--- WAIC computation
Waic <- function(stanfit) {
  loglik <- extract (stanfit, "logLik")$logLik
  n1 <- nrow(loglik)  # number of samples
  n2 <- ncol(loglik)  # number of data points
  vars <- colSums((loglik - matrix(colMeans(loglik), n1, n2, byrow = T))^2) / (n1 - 1)
  pwaic <- sum(vars)  # effective parameter number
  lpd <- sum(log(colMeans(exp(loglik))))  # log pointwise predictive density
  waic <- -2 * (lpd - pwaic)
  return(list(WAIC = waic, Pwaic = pwaic, LPD = lpd))
}

# #--- debug
# model <- stan_model(file = "~/Dropbox/ADSP/bayes/hir_null.stan", 
# 	     model_name='GLM (Hierarchical)', verbose = FALSE)
# 
# dataList1 <- c("N", "K", "L", "Ad", "Age", "Sex", "Family")
# Snp <- rep(0, N)
# dataList2 <- c(dataList1, "Snp")
# fit <- sampling(model, data = dataList2, iter = 1000, chains = 3, refresh = 1e2)
# 
# print(fit, pars = c("pAge", "pSex"), probs = c(0.1, 0.5, 0.9))

#--- STAN model
hir_model_null <- stan_model(file = "~/Dropbox/ADSP/bayes/hir_null.stan", 
         	    model_name='GLM (Hierarchical)', verbose = FALSE)
hir_model <- stan_model(file = "~/Dropbox/ADSP/bayes/hir.stan", 
         	    model_name='GLM (Hierarchical)', verbose = FALSE)

#--- null model
nchains <- 4

dataList1 <- c("N", "K", "L", "Ad", "Age", "Sex", "Family")
fit1 <- sampling(model1, data=dataList1, iter = 1, chains = 1)  # init
fit2 <- mclapply(1:nchains, mc.cores = nchains, FUN = function (chain) {
  stan(fit = fit1, data = dataList1, iter = 2000, 
       chains = 1, seed = 0, chain_id = chain, refresh = 0)})
fit <- sflist2stanfit(fit2)

waic1 <- Waic(fit)$WAIC

#--- full model
load("~/Dropbox/ADSP/R/vdata_21loci.rdt")
n.snp <- nrow(geno)
waic2 <- rep(0, nrow(geno))
names(waic2) <- rownames(geno)

cat(date(), "--- start \n")
#---------------------------------------------------------------------
for (i in 1:nrow(geno)) {
  cat(i, "in", nrow(geno), rownames(geno)[i], "\n")
  Snp = geno[i, ]
  dataList2 <- c(dataList1, "Snp")

  fit1 <- sampling(model2, data=dataList2, iter = 1, chains = 1)  # init
  fit2 <- mclapply(1:nchains, mc.cores = nchains, FUN = function (chain) {
    stan(fit = fit1, data = dataList2, iter = 2000, 
         chains = 1, seed = i, chain_id = chain, refresh = 0)})
  fit <- sflist2stanfit(fit2)

  waic2[i] <- Waic(fit)$WAIC
}
save(waic1, waic2, file = "~/Dropbox/ADSP/R/waic_21loci.rdt")
#---------------------------------------------------------------------
cat(date(), "--- finish \n")

