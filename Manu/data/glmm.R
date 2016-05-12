library(dplyr)
library(ggplot2)
library(rstan)
library(QTLRel)

# Why is GLMM by Stan better than LMM by QTLRel?

# 
rm(list = ls())
hpc <- "/data/xwang/Adsp"
github <- "~/Dropbox/GitHub/Adsp"

# 
glmm = list() # output

setwd(github)
save(glmm, file = "./data/glmm.rdt")

#
setwd(hpc)
setwd("~/Desktop/")

load("./GWAS/optimizing.rdt")
glm <- do.call(rbind, glm.fit)
rm(glm.fit)

load("./GWAS/QtlRel.rdt") 
lmm <- do.call(rbind, lmm.fit[1:22])
rm(lmm.fit)

all(glm$UID == lmm$UID)

dLOD = lmm$LOD - glm$LOD
summary(dLOD)

(index = which(dLOD > 7))

fit = cbind(glm[index, ], lmm[index, 17:24])
rownames(fit) = NULL

var = fit$UID

setwd(hpc)

geno <- lapply(1:22, function(chr) {
  cat(chr, "\n")
  load(paste0("./golden/geno_chr", chr, ".rdt"))
  geno[rownames(geno) %in% var, ]
})

geno <- do.call(rbind, geno)

#
setwd(github)
load("./data/glmm.rdt")
for(obj in names(glmm)) assign(obj, glmm[[obj]])

colnames(fit)[21] = "pSNP"
colnames(fit)[23] = "GLMM"
colnames(fit)[30] = "LMM"

fit$dLOD = fit$LMM - fit$GLMM

fit = fit[fit$dLOD > 8, ]

table(fit$GLMM > 15)

fit = fit[order(fit$dLOD), ]

qplot(MAF, dLOD, data = fit, size = GLMM) # graph
qplot(MAF, dLOD, data = fit, size = abs(pSNP))
qplot(GLMM, dLOD, data = fit, size = abs(pSNP))

fit$Sign = "Risky"
fit$Sign[fit$pSNP < 0] = "Protective"

pdf("./Manu/glmm1.pdf", width = 6, height = 4, family = "Helvetica")

ggplot(fit, aes(x = MAF, y = dLOD)) + 
  geom_point(aes(size = GLMM, colour = Sign), shape = 111) + 
  scale_size(range = c(2, 15)) +
  theme_bw() + xlab("MAF") + ylab("d-LOD") + 
  scale_color_manual(values = c("dodgerblue3", "firebrick1")) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(size = .5),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key = element_blank()) 

dev.off()

(t1 = fit[which.max(fit$dLOD), ])
(t2 = fit[which.max(fit$GLMM), ])

xx = mutate(mdata, geno = geno[t1$UID, ])
xx = mutate(mdata, geno = geno[t2$UID, ])

qplot(x = AD2, data = xx, geom = "bar", fill = as.factor(geno), position = "fill")
table(xx$geno, xx$AD2)

mycol <- c("dodgerblue3", "firebrick1")

pdf("./Manu/glmm2.pdf", width = 6, height = 2.5, family = "Helvetica")

ggplot(xx, aes(x = AD2, fill = as.factor(geno))) + 
  geom_bar(position = "fill", width = 0.65) +
  theme_bw() + xlab("") + ylab("Frequency")  + coord_flip() +
  scale_fill_manual(values = mycol, 
                    breaks=c("0", "1", "2"), 
                    labels=c("0/0", "0/1", "1/1")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'grey30'),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.key = element_blank()) 

dev.off()

vep[vep$UID == t1$UID, ]

table(mdata$AD2, geno[t1$UID, ])

for (i in 1:117) {
  if (fit$UID[i] %in% rownames(geno)) {
    print(fit$UID[i])
    xx = mutate(mdata, geno = geno[fit$UID[i], ])
    print(table(xx$geno, xx$AD2))
  }
}

# Stan's cutpoints estimation do not change much for the 117

load("./data/sampling.rdt")
glm_stan <- stan_model("./GLM/glm(151024).stan")

cov <- data$mdata[c("Age", "Sex", "Apoe2", "Apoe4")]
cov$Sex <- Sex <- as.integer(cov$Sex) - 1
Ad <- as.numeric(data$mdata$AD2)
init <- data$init

cov$Var <- geno[t1$UID, ]
Sigma <- data$kinship[[paste0("no_chr", 19)]]

cov$Var <- geno[t2$UID, ]
Sigma <- data$kinship[[paste0("no_chr", 11)]]

Sigma[Sigma < 0] <- 0
L <- t(chol(Sigma))

dat = list(N = 576, K = 4, D = 5, cov = cov, L = L, Ad = Ad)

fit = optimizing(glm_stan, verbose = FALSE, algorithm = "LBFGS", init = init, data = dat)
loglik = sum(fit$par[paste0("lp[", 1:576, "]")])
fit$par[c("c[1]", "c[2]", "c[3]")]

geno = geno[! rownames(geno) == "", ]

cuts <- lapply(1:nrow(geno), function(x) {
  dat0 = dat
  dat0$cov$Var = geno[x, ]
  
  chr = gsub("-.*", "", rownames(geno)[x])
  Sigma <- data$kinship[[paste0("no_chr", chr)]]
  Sigma[Sigma < 0] <- 0
  L <- t(chol(Sigma))
  dat0$L = L
  
  fit = optimizing(glm_stan, verbose = FALSE, algorithm = "LBFGS", init = init, data = dat0)
  fit$par[c("c[1]", "c[2]", "c[3]")]
})

cuts = do.call(rbind, cuts)

boxplot(cuts)

which(rownames(geno) == t1$UID)
cuts[100, ]

# QTLRel: how much LOD change by another numerical coding?

Y <- mdata$AD1

Y[Y == 0.25] = 0.33
Y[Y == 0.50] = 0.66

Y[Y == 0.25] = 0.20
Y[Y == 0.50] = 0.60

X <- mdata[c("Age", "Sex", "Apoe2", "Apoe4")]

KS <- data$kinship$autosome
KS <- data$kinship[[paste0("no_chr", 19)]]
KS[KS < 0] <- 0

o <- estVC(y=Y, x=X, v=list(AA=KS, DD=NULL, HH=NULL, AD=NULL, MH=NULL, EE=diag(576)))

g1 <- geno[t1$UID, ] + 1 # QTLRel use 1, 2, 3
scanOne(Y, x = X, gdat = as.matrix(g1), numGeno = TRUE, vc = o)
