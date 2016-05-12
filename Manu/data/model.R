library(rstan)
library(reshape)
library(dplyr)

rm(list = ls())
setwd("~/Dropbox/GitHub/Adsp")
load("./Manu/model.rdt") 

names(model)

beta = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]")
name = c("Age", "Sex", "APOE/e2", "APOE/e4", "Var")

fit = model$sampling
sample <- as.data.frame(fit)[beta]
print(fit, pars = beta)

mode <- apply(sample, 2, function(z) { dens <- density(z); dens$x[which.max(dens$y)] })
mean = colMeans(sample)
se1 = apply(sample, 2, function(xx) quantile(xx, 0.05))
se2 = apply(sample, 2, function(xx) quantile(xx, 0.95))
  
dt = data.frame(feature = name, mode, mean, se1, se2)
dt$feature = factor(dt$feature, levels = name)

col4 <- c("dodgerblue3", "darkorchid2", "chartreuse3", "firebrick1")
col5 <- c("darkorchid2", "dodgerblue3", "gold1", "chartreuse3", "firebrick1")
col6 <- c("grey70", "firebrick1", "dodgerblue3", "gold1", "chartreuse3", "darkorchid2")

pdf("./Manu/covar.pdf", width = 5, height = 3.5, family = "Helvetica")

ggplot(dt[1:4, ], aes(x = feature, y = mode)) + 
  geom_errorbar(aes(ymax = se1, ymin = se2, color = feature), size = 1, width = .3) + 
# geom_point(aes(size = abs(mode), color = feature)) + 
  geom_point(aes(color = feature), size = 7) + 
  scale_color_manual(values = col4) + coord_flip() +
  geom_hline(y = 0, linetype = 2, size = 1) +
  theme_bw() + xlab("") + ylab("Effect Size") + ylim(-1.5, 1.5) +
  theme(panel.border = element_rect(size = 1, color = "grey30"),
        axis.text = element_text(size = 13),
      	axis.title = element_text(size = 15),
        legend.position = "none")

dev.off()

models <- c("age_sex", "age_apoe2", "age_apoe4", "sex_apoe2", "sex_apoe4", "apoe2_apoe4")

df = lapply(models, function(idx) {
  fit = model[[idx]]
  sample <- as.data.frame(fit)[beta]
  
  mode <- apply(sample, 2, function(z) { dens <- density(z); dens$x[which.max(dens$y)] })
  se1 = apply(sample, 2, function(xx) quantile(xx, 0.05))
  se2 = apply(sample, 2, function(xx) quantile(xx, 0.95))
  
  dt = data.frame(feature = name, mode, se1, se2)
  dt$feature = factor(dt$feature, levels = name)
  dt$model = idx
  dt
})

dt = do.call(rbind, df)

dt$feature = gsub("Var", "Interaction", dt$feature)
dt$feature = factor(dt$feature, levels = c("Age", "Sex", "APOE/e2", "APOE/e4", "Interaction"))

mylab <- c("Age by Sex", "Age by APOE/e2", "Age by APOE/e4", "Sex by APOE/e2", "Sex by APOE/e4", "APOE/e2 by APOE/e4")
dt$model <- factor(dt$model, levels = models, labels = mylab)

names(dt) = gsub("feature", "Feature", names(dt))
names(dt) = gsub("mode", "Mode", names(dt))

pdf("./Manu/covar_inter.pdf", width = 10, height = 5, family = "Helvetica")

ggplot(dt, aes(x = Feature, y = Mode)) + 
  geom_errorbar(aes(ymax = se1, ymin = se2, color = Feature), width = 0.2) + 
  geom_point(aes(size = abs(Mode), color = Feature)) + facet_grid(. ~ Model) +
  geom_hline(y = 0, linetype = 2, size = .3) + coord_flip() +
  theme_bw() + xlab("") + ylab("") + ylim(-2.5, 2.0) +
  scale_color_manual(values = col5) + 
  theme(panel.border = element_rect(size = 1, color = "grey30"),
        axis.text = element_text(size = 13),
      	axis.title = element_text(size = 15, vjust = 1),
        legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 12),
        legend.key = element_blank())

dev.off()
  
# --- QTLRel

KS <- kinship$autosome
vc <- estVC(y = data0[, "AD1"], x = data0[, c("Age", "Sex", "APOE")], 
        v = list(AA = KS, EE = diag(nrow(data0)), DD = NULL, HH = NULL, AD = NULL, MH = NULL))

snp <- sample(1:3, nrow(data0), replace = T)
qtlrel0 <- scanOne(y = data0[, "AD1"], x = data0[, c("Age", "Sex", "APOE")], 
                   vc = vc, gdat = as.matrix(snp))
qtlrel0$p
qtlrel0$v
qtlrel0$parameters

eff.lm0 <- summary(lm0)$coefficients[-1, ]
eff.lmer0 <- summary(lmer0)$coefficients[-1, ]
eff.clm0 <- summary(clm0)$coefficients[-c(1:3), ]
eff.clmm0 <- summary(clmm0)$coefficients[-c(1:3), ]

features <- c("Age", "Sex", "APOE22", "APOE23", "APOE24", "APOE34")
eff.clm <- data.frame(feature = factor(rep(features, 2), levels = features), 
                      model = factor(rep(c("CLM", "CLMM"), each = 6)), 
                      effect = c(eff.clm0[, 1], eff.clmm0[, 1]),
                      se = c(eff.clm0[, 2], eff.clmm0[, 2]))

ggplot(eff.clm, aes(x = as.numeric(model), y = effect, color = model)) +  
  geom_point(aes(shape = model), size = 5) + 
  geom_errorbar(aes(ymax = effect + se, ymin = effect - se), width = 0.2) + 
  facet_grid( . ~ feature) +
  theme_bw() + xlab("") + ylab("Effect") + 
  scale_color_manual(values = c("firebrick1", "dodgerblue3")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, face = "bold", vjust = 2),
        legend.title = element_blank(), 
        legend.text = element_text(size = 10, face = "bold"),
        legend.key = element_blank()) 
