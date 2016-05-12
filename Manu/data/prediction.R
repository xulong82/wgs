library(xlsx)
library(dplyr)
library(ggplot2)
library(stargazer)
library(MASS)

rm(list = ls())
setwd("~/Dropbox/GitHub/Adsp")

options(stringsAsFactors = F)

load("./data/mdata.rdt")
load("./data/glmList.rdt")
load("./Manu/data.rdt")

gwas <- data$top$gwas
gwas <- filter(gwas, LOD > 15 & CHR %in% c(1:22))
gwas <- filter(gwas, UID %in% PEAK)

geno <- glmList$geno

table(gwas$UID %in% rownames(geno))
uid <- intersect(gwas$UID, rownames(geno))

geno <- geno[uid, ]
gwas <- gwas[match(uid, gwas$UID), ]

all(colnames(geno) == mdata$SRR)
effectCum <- sweep(geno, 1, gwas$pSnp, "*")
mdata$Cum <- colSums(effectCum)

meta_extra <- read.xlsx("./docs/adsp_subject_and_experiment_metadata.xlsx", startRow = 2, sheetIndex = 1)
mdata$Center <- meta_extra$LSAC[match(mdata$SRR, meta_extra$SRR)]
mdata$Ethnicity <- meta_extra$Ethnicity[match(mdata$SRR, meta_extra$SRR)]
mdata$Ethnicity[mdata$Ethnicity == "0"] = "Not-Latino"
mdata$Ethnicity[mdata$Ethnicity == "1"] = "Latino"
mdata$Ethnicity[mdata$Ethnicity == "null"] = "Null"

mdata$Race[mdata$Race == "4"] = "Black"
mdata$Race[mdata$Race == "5"] = "White"
mdata$Race[mdata$Race == "6"] = "Other"
mdata$Race[mdata$Race == "null"] = "Null"

mdata$Race = factor(mdata$Race, levels = c("Black", "White", "Other", "Null"))

table(mdata$AD2, mdata$Center) %>% as.data.frame.matrix
stargazer(table(mdata$AD2, mdata$Center) %>% as.data.frame.matrix, summary = F)

summary(gwas$pSnp)
summary(mdata$Cum)

geno_def = geno[, mdata$AD2 == "Definite"]
geno_indef = geno[, mdata$AD2 != "Definite"]
which(rowSums(geno_def) & ! rowSums(geno_indef)) # no one only in definite category

cov = mdata[c("Age", "Sex", "Apoe2", "Apoe4")]

step.dat = data.frame(y = mdata$AD1, t(geno))
step.dat = data.frame(y = mdata$AD1, t(effectCum))
step.dat = data.frame(y = mdata$AD1, t(effectCum), cov)

null = lm(y ~ 1, data = step.dat)
full = lm(y ~ ., data = step.dat)

step.fit = step(null, scope = list(upper=full), data = step.dat, direction="both")

AIC = step.fit$anova$AIC
names(AIC) = step.fit$anova$Step

plot(step.fit$anova$AIC)

coef = step.fit$coefficients

step.var = names(coef)[grep("^X", names(coef))]
step.var = gsub("\\.", "-", step.var)
step.var = gsub("X", "", step.var)

mdata$Cum <- colSums(geno[step.var, ])
mdata$Cum <- colSums(effectCum[step.var, ])

cor(mdata$Cum, mdata$AD1) # additive? interactive?

mycol <- c("grey70", "dodgerblue3", "chartreuse3", "firebrick1")

pdf("temp/cum(15).pdf", width = 7, height = 5)

ggplot(mdata, aes(x = AD2, y = Cum, fill = Race)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(shape = Ethnicity), position = position_jitter(width = 0.3)) +
  theme_bw() + xlab("") + ylab("") +
  scale_fill_manual(values = mycol) +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'grey30'),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key = element_blank()) 

dev.off()

# what are the variants combinations that only appear in definite category?

filter(mdata, Add > 100) # these combination all definite disease
SRR <- filter(mdata, Add > 100)$SRR # highly risky individuals
genoSRR <- geno_select[, SRR]
genoSRR <- genoSRR[rowSums(genoSRR) != 0, ]

# Use this variants panel, calculate the score to predict the risk

rowSums(genoSRR & 1) %>% summary
colSums(genoSRR & 1) %>% summary # they have 67 to 127 variants

profiles = apply((genoSRR & 1) + 0, 2, function(x) paste(x, collapse = ""))
table(profiles) # they each have unique variants combination 
