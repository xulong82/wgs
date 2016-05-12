library(arm)
library(lme4)
library(rstan)
library(ordinal)
library(emma)
library(QTLRel)
library(lmtest)
library(ggplot2)
library(dplyr)

rm(list = ls())
hpc <- "/data/xwang/ADSP"
github <- "~/Dropbox/GitHub/Adsp"

setwd(github)

load("./data/kinship.rdt")

load("./data/mdata.rdt")
load("./diagnosis/chr1pList.rdt")

glm1 <- stan_model(file = "./GLM/glm.stan", model_name='GLMM')  # Factorial APOE
glm2 <- stan_model(file = "./GLM/glm_new.stan", model_name='GLMM')  # Apoe3/2/4

dat1 <- mdata[c("AD2", "Age", "Sex", "APOE")] %>% as.list
dat1 <- within(dat1, { Sex <- as.numeric(Sex); Ad <- as.numeric(AD2); Apoe = as.numeric(APOE) })
dat1$AD2 <- dat1$APOE <- NULL
dat1 <- within(dat1, { N = 576; K = 4; KS = KS; Snp = chr1p$one$geno })
fit1 <- optimizing(glm1, algorithm = "LBFGS", data = dat1)
fit1$par[c("pAge", "pSex", "pSnp", "pApoe22", "pApoe23", "pApoe24", "pApoe34")] 
sum(fit1$par[paste("logLik[", 1:576, "]", sep = "")])

dat2 <- mdata[c("AD2", "Age", "Sex", "Apoe2", "Apoe4")] %>% as.list
dat2 <- within(dat2, { Sex <- as.numeric(Sex); Ad <- as.numeric(AD2); rm(AD2) })
dat2 <- within(dat2, { N = 576; K = 4; KS = KS; Snp = chr1p$one$geno })

fit2 <- optimizing(glm2, algorithm = "LBFGS", data = dat2)
fit2$par[c("pAge", "pSex", "pSnp", "pApoe2", "pApoe4")] 
sum(fit2$par[paste("logLik[", 1:576, "]", sep = "")])

mcmc1 <- sampling(glm1, data = dat1, warmup = 3e2, iter = 6e2, chains = 3)  # why so long???
mcmc1 <- as.data.frame(mcmc1)[c("pAge", "pSex", "pApoe2", "pApoe4")]
mode <- apply(mcmc1, 2, function(z) { dens <- density(z); dens$x[which.max(dens$y)] })
ci95 <- summary(mcmc1)$summary[c("pAge", "pSex", "pApoe2", "pApoe4"), c("2.5%", "97.5%")]

mdata$Snp <- chr1p$one$geno
(lm <- lm(formula = AD1 ~ Age + Sex + APOE + Snp, mdata))
(lm <- lm(formula = AD1 ~ Age + Sex + Apoe2 + Apoe4 + Snp, mdata))

(lmer <- lmer(formula = AD1 ~ Age + Sex + APOE + Snp + (1 | Family.ID), mdata, REML = F))
(clm <- clm(formula = AD2 ~ Age + Sex + APOE + Snp, data = mdata))
(clmm <- clmm(formula = AD2 ~ Age + Sex + APOE + Snp + (1 | Family.ID), mdata))

X <- data0[, "AD1"]
Y <- data0[, c("Age", "Sex", "APOE")]
KS <- kinship$autosome

vc <- estVC(y = Y, x = X, v = list(AA=KS, EE=diag(576), DD=NULL, HH=NULL, AD=NULL, MH=NULL))
qtlrel <- scanOne(y = Y, x = X, vc = vc, gdat = as.matrix(snp))
qtlrel0$p
qtlrel0$v
qtlrel0$parameters

setwc(hpc) # Apoe locus - 200k
load("biVar/biGeno_chr19.rdt")
load("Meta/meta_chr19.rdt")
meta = meta[meta$MAF > 0.01 & meta$HET < 0.99, ]
meta = meta[meta$POS > 45330000 & meta$POS < 45530000, ]
geno = biGeno[meta$UID, ]
Apoe = list(geno = geno, meta = meta)
save(Apoe, file = paste0(github, "/diagnosis/ApoeList.rdt"))

setwd(hpc) # rs139401485 locus - 200k
load("biVar/biGeno_chr1.rdt")
load("Meta/meta_chr1.rdt")
meta = meta[meta$MAF > 0.01 & meta$HET < 0.99, ]
meta = meta[meta$POS > 1570415 & meta$POS < 1770415, ]
geno = biGeno[meta$UID, ]
chr1p = list(geno = geno, meta = meta)
save(chr1p, file = paste0(github, "/diagnosis/chr1pList.rdt"))

mdata$geno <- as.factor(chr1p$one$geno)
mycol <- c("grey70", "dodgerblue3", "firebrick1")
pdf("./diagnosis/geno.pdf", width = 4.5, height = 3, family = "Helvetica")
ggplot(mdata, aes(x = AD2, fill = geno)) + 
  geom_bar(position = "fill", width = 0.65) +
  theme_bw() + xlab("") + ylab("")  + coord_flip() +
  scale_fill_manual(values = mycol, breaks=c("0", "1", "2"), labels=c("0/0", "0/1", "1/1")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'grey30'),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.key = element_blank()) 
dev.off()
