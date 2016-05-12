library(lme4)
library(ordinal)
library(lmtest)
library(parallel)
library(dplyr)

rm(list = ls())
hpc <- "/data/xwang/ADSP"
github <- "~/Dropbox/GitHub/Adsp"

arg <- commandArgs(TRUE) 
chr <- paste("chr", arg, sep = "")

setwd(github)
load("./data/mdata.rdt")

setwd(hpc)
load(paste("biVar/biGeno_", chr, ".rdt", sep = "")) 
load(paste("Meta/meta_", chr, ".rdt", sep = "")) 

uid <- meta$UID[meta$MAF > 0.01 & meta$HET < 0.99]
geno <- biGeno[uid, ]

n.core <- 20 
u.core <- round(nrow(geno) / n.core)
idx1 <- ((1:n.core) - 1) * u.core + 1 
idx2 <- c((1:(n.core - 1)) * u.core, nrow(geno))
gList <- mclapply(1:20, function(x) geno[idx1[x]:idx2[x], ], mc.preschedule = F)

dat0 <- mdata[c("AD1", "AD2", "Sex", "Age", "Apoe2", "Apoe4", "Family.ID")]
lmer0 <- lmer(AD1 ~ Sex + Age + Apoe2 + Apoe4 + (1 | Family.ID), dat0, REML = FALSE) %>% logLik
clmm0 <- clmm(AD2 ~ Sex + Age + Apoe2 + Apoe4 + (1 | Family.ID), dat0) %>% logLik

myLmm <- function (x) {
  fit <- matrix(nrow = nrow(x), ncol = 4, 
    dimnames = list(rownames(x), c("LMER_SIZE", "LMER_LOD", "CLMM_SIZE", "CLMM_LOD")))

  for (i in 1:nrow(x)) {
    dat1 <- within(dat0, VAR <- x[i, ])
    lmer <- lmer(AD1 ~ Sex + Age + Apoe2 + Apoe4 + VAR + (1 | Family.ID), dat1, REML = FALSE)
    clmm <- clmm(AD2 ~ Sex + Age + Apoe2 + Apoe4 + VAR + (1 | Family.ID), dat1)

    lmer.coef <- summary(lmer)$coefficients
    clmm.coef <- summary(clmm)$coefficients
    lmer.size <- ifelse("VAR" %in% rownames(lmer.coef), lmer.coef["VAR", "Estimate"], 0)
    clmm.size <- ifelse("VAR" %in% rownames(clmm.coef), clmm.coef["VAR", "Estimate"], 0)
    fit[i, ] <- c(lmer.size, logLik(lmer), clmm.size, logLik(clmm))
  }

  fit[, "LMER_LOD"] <- 2 * (fit[, "LMER_LOD"] - lmer0)
  fit[, "CLMM_LOD"] <- 2 * (fit[, "CLMM_LOD"] - clmm0)

  return(fit)
}

y <- do.call(rbind, mclapply(gList, myLmm, mc.cores = 20))
assign(chr, y)
save(list = chr, file = paste0("LMM/", chr, ".rdt"))

