library(ggplot2)

rm(list = ls())
setwd("~/Dropbox/GitHub/Adsp")
load("./Manu/sampling.rdt")

opt = sampling$opt
waic = sampling$waic
summary = sampling$summary

opt$waic = max(waic) - waic
opt$POS = opt$POS * 1e-3

opt$opt_mean = opt$"beta[5]"
opt$spl_mean = sapply(summary, function(x) x["beta[5]", "mean"])

opt$SNP = factor(opt$SNP, levels = c("Y", "N"))

mycol = rep("chartreuse3", nrow(opt))
mycol[opt$SNP == "N"] = "firebrick1"

pdf("./Manu/scatter.pdf", width = 6, height = 4)

ggplot(opt, aes(x = opt_mean, y = spl_mean)) + 
  geom_point(aes(colour = SNP, size = MAF)) +
  geom_abline(intercept = 0, slope = 1) + 
  scale_color_manual(values = c("dodgerblue3", "firebrick1")) + 
  theme_bw() + xlab("Effect - Optimizing") + ylab("Effect - Sampling") +
  theme(panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key = element_blank()) 

dev.off()

plot(x = opt$"beta[5]", y = mc.mean, type = "p", pch = 23, cex = 1.0, bg = mycol)

pdf(file = "./Manu/waic.pdf", width = 7, height = 10)

close.screen(all.screens = TRUE)
split.screen(rbind(c(0.1, 0.9, 0.1, 0.5), c(0.1, 0.9, 0.5, 0.9)))

screen(1)
par(mar = c(0, 0, 0, 0))

plot(x = opt$POS, y = opt$LOD, type = "p", pch = 23, cex = 1.0, bg = mycol, 
     main = "", xlab = "", ylab = "", ylim = c(0, 25), axes = F)

axis(1, at = c(min(opt$POS), max(opt$POS)), las = 1, lwd = 2) 
mtext(paste("Chromosome 1", "Position (Kb)", sep=" "), side = 1, line = 2.5, font = 2)

axis(2, at = c(0, 10, 20), labels = c(0, 10, 20), las = 1, lwd = 2) 
mtext("LOD", side = 2, at = 10, line = 2, font = 2)

hit = opt[which.max(opt$LOD), ]
points(hit$POS, hit$LOD, pch = 5, cex = 2.0, lwd = 2.0, col = "blue")
text(hit$POS, hit$LOD, labels = hit$ID, pos = 3, offset = 1, font = 2)

screen(2)
par(mar = c(0, 0, 0, 0))

plot(x = opt$POS, y = opt$waic, type = "p", pch = 23, cex = 1.0, bg = mycol, 
     main = "", xlab = "", ylab = "", ylim = c(0, 30), axes = F)

# abline(0, 0, lwd = 2, col = "black")

axis(2, at = c(0, 10, 20, 30), labels = c(0, 10, 20, 30), las = 1, lwd = 2) 
mtext("WAIC", side = 2, at = 15, line = 2, font = 2)

points(hit$POS, hit$waic, pch = 5, cex = 2.0, lwd = 2.0, col = "blue")
text(hit$POS, hit$waic, labels = hit$ID, pos = 3, offset = 1, font = 2)

close.screen(all.screens = TRUE)
dev.off()
