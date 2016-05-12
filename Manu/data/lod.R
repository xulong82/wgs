library(ggplot2)
library(pheatmap)

rm(list = ls())
setwd("~/Dropbox/GitHub/Adsp")
load("./Manu/lod.rdt")

lod = lod$lod2

lapply(lod, cor)

(mycor = sapply(lapply(lods, cor), function(x) x[1, -1]))
colnames(mycor) = c("5%", "1%", "0.5%", "0.1%", "0.05%", "0.01%", "0.005%")
mycor = mycor[, 7:1]

dt = melt(mycor)
names(dt) = c("Models", "Quantile", "value")

pdf("./Manu/lod-cor.pdf", width = 6, height = 4)

ggplot(dt, aes(x = Quantile, y = value, group = Models)) + 
  geom_line(aes(color = Models), size = 1) + geom_point(size = 3) +
  theme_bw() + xlab("") + ylab("Correlation coefficients") +
  scale_color_manual(values = c("dodgerblue3", "firebrick1", "chartreuse3")) + 
  theme(panel.border = element_rect(size = 1, color = "grey30"),
        axis.text = element_text(size = 13, angle = 90),
        axis.title = element_text(size = 16),
        legend.key = element_blank())

dev.off()

pdf("./Manu/lod-heatmap.pdf", width = 2, height = 2)

pheatmap(cor(lods[[1]]), display_number = T, treeheight_col = 0, treeheight_row = 0, legend = F)

dev.off()

setwd("/data/xwang/Adsp")

load("./GWAS/lodALL.rdt") 

lodALL = lod[c("GLMM", "LMM")]
lodALL = lod[c("GLMM", "LMM", "LMER", "CLMM")]

x = lodALL$LMM - lodALL$GLMM

setwd("~/Dropbox/GitHub/Adsp")

png("./Manu/pairs.png", width = 2000, height = 2000, res = 300)

pairs(lodALL, cex = .1)

dev.off()

png("./Manu/glm-lmm.png", width = 1200, height = 1000, res = 200)

ggplot(lodALL, aes(x = GLMM, y = LMM)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw() + xlab("LOD - GLMM") + ylab("LOD - LMM") +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "firebrick1") +
  xlim(c(0, 36)) + ylim(c(0, 36)) +
  theme(panel.border = element_rect(size = 1, color = "grey30"),
  axis.text = element_text(size = 13),
  axis.title = element_text(size = 16))

dev.off()

