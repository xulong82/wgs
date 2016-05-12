#  Bayesian inference on ADSP project
#  -------------------------------------------------------------------
#  Random effect: Genotype-based Relatedness

library(rstan)
library(parallel)

rm(list = ls())

#--- META data
load("~/Dropbox/ADSP/R/mdata.rdt")
load("~/Dropbox/ADSP/R/vdata_21loci.rdt")
data1 <- mdata[, c("AD1", "Sex", "Age")] # Sex: 0-M, 1-F
data1$AD2 <- factor(data1$AD1, levels = c("0", "0.25", "0.5", "1"), ordered = T) 

N <- nrow(data1)  # sample number
K <- length(unique(data1$AD1))  # AD categories
Ad <- as.numeric(data1$AD2)
Age <- data1$Age
Sex <- as.numeric(data1$Sex)

load("~/Dropbox/ADSP/R/kinship2.rdt")
KS <- kin2

#--- STAN model: debug
model <- stan_model(file = "~/Dropbox/ADSP/Bayes/kinship_null.stan", 
                    model_name='ordered_logitistic model (Null)', verbose = FALSE)
dataList <- c("N", "K", "KS", "Ad", "Age", "Sex")
fit <- sampling(model, data = dataList, iter = 1000, chains = 1, refresh = 10)

#--- STAN model
glm_model_null <- stan_model(file = "~/Dropbox/ADSP/bayes/glm_null.stan", 
                     model_name='ordered_logitistic model (Null)', verbose = FALSE)
glm_model <- stan_model(file = "~/Dropbox/ADSP/bayes/glm.stan", 
                     model_name='ordered_logitistic model (Full)', verbose = FALSE)

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

#--- Null model
nchains <- 4

dataList1 <- c("N", "K", "KS", "Ad", "Age", "Sex")
fit1 <- sampling(model1, data=dataList1, iter = 1, chains = 1)  # init
fit2 <- mclapply(1:nchains, mc.cores = nchains, FUN = function (chain) {
  stan(fit = fit1, data = dataList1, iter = 2000, 
       chains = 1, seed = 0, chain_id = chain, refresh = 0)})
fit <- sflist2stanfit(fit2)

waic1 <- Waic(fit)$WAIC

#--- Full model
# loci <- read.delim("~/Dropbox/ADSP/21loci.txt", stringsAsFactors = F)
# loci <- loci[loci$SNP %in% rownames(geno), ]
# geno <- geno[loci$SNP, ]
# save(geno, file = "~/Dropbox/ADSP/R/vdata_21loci.rdt")

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
save(waic1, waic2, file = "~/Dropbox/ADSP/bayes/waic_21loci.rdt")
#---------------------------------------------------------------------
cat(date(), "--- finish \n")
