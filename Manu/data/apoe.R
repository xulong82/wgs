library(dplyr)
library(ggplot2)
library(stargazer)

rm(list = ls()) 
setwd("~/Dropbox/GitHub/Adsp")

load("./data/hg19.rdt")
load("./data/apoe.rdt")

# APOE allele frequencies

(47 + 4) / 576 / 2
(243 + 650) / 576 / 2
(208) / 576 / 2

chr <- apoe$hg19$CHR
pos <- apoe$hg19$POS
  
min.pos <- pos - 250e3
max.pos <- pos + 250e3 
  
hg19 = filter(hg19, CHR == chr, POS > min.pos, POS < max.pos)
hg19 = hg19[order(hg19$POS), ]

hg19 = hg19[hg19$NAME != "APOC4-APOC2", ]

null = apoe$opt$null 
full = apoe$opt$full

all(rownames(null) == full$UID)
full$LOD2 = null$LOD
full$pVar = null$"beta[3]"

pos <- pos * 1e-3
min.pos <- min.pos * 1e-3 
max.pos <- max.pos * 1e-3

full$POS <- full$POS * 1e-3

hg19$POS <- hg19$POS * 1e-3
hg19$START <- hg19$START * 1e-3
hg19$END <- hg19$END * 1e-3
  
hit <- full[full$LOD2 == max(full$LOD2), ]

pdf(file = "./Manu/apoe.pdf", width = 7, height = 10)

close.screen(all.screens = TRUE)

split.screen(rbind(c(0.1, 0.9, 0.1, 0.3), c(0.1, 0.9, 0.3, 0.6), c(0.1, 0.9, 0.6, 0.9)))

screen(1)
par(mar = c(0, 0, 0, 0))

plot(x = full$POS, y = rep(0, length(full$POS)), ylim = c(-12, 0), type = "n", axes = F)
# box()

for (i in 1:nrow(hg19)) {  # plot the genes
  adj.arrow <- -(i %% 5 + 1) * 2 
  adj.text <- -(i %% 5 + 1) * 2 - 1 
  
  text(hg19[i, ]$POS, adj.text, labels = hg19[i, ]$NAME, cex = .7, pos = 1, offset = 0, font = 2)
  arrows(max(hg19[i, ]$START, min.pos), adj.arrow, min(hg19[i, ]$END, max.pos), adj.arrow, 
             length = 0.05, lwd = 2, code = 2, lty = "solid", col = "darkgreen")
}

axis(1, at = c(min.pos, pos, max.pos), labels = round(c(min.pos, pos, max.pos)), las = 1, lwd = 2) 
mtext(paste("Chromosome", chr, "Position (Kb)", sep=" "), side = 1, line = 2.5, font = 2)

screen(2)
par(mar = c(0, 0, 0, 0))

mycol = rep("chartreuse3", nrow(full))
mycol[full$pVar > 0] = "firebrick1"

plot(x = full$POS, y = null$LOD, type = "p", pch = 23, cex = 1.0, bg = mycol,
     main = "", xlab = "", ylab = "", ylim = c(0, 25), axes = F)

axis(2, at = c(0, 10, 20), labels = c(0, 10, 20), las = 1, lwd = 2) 
# mtext("LOD", side = 2, at = 10, line = 2, font = 2)

hit = full[full$ID %in% c("rs429358", "rs7412"), ]
points(hit$POS, hit$LOD2, pch = 5, cex = 2.0, lwd = 2.0, col = "blue")
# points(hit$POS, hit$LOD2, pch = 5, cex = 2.0, lwd = 2.0, col = "blue")
# text(hit$POS, hit$LOD2, labels = hit$ID, pos = 3, offset = 1, font = 2)

lines(c(min.pos, max.pos), c(15, 15), lty = "dotted", lwd = 1, col = "blue")

screen(3)
par(mar = c(0, 0, 0, 0))

mycol = rep("chartreuse3", nrow(full))
mycol[full$"beta[5]" > 0] = "firebrick1"

plot(x = full$POS, y = full$LOD, type = "p", pch = 23, cex = 1.0, bg = mycol, 
     main = "", xlab = "", ylab = "", ylim = c(0, 25), axes = F)

hit = full[full$ID %in% c("rs429358", "rs7412"), ]
points(hit$POS, hit$LOD, pch = 5, cex = 2.0, lwd = 2.0, col = "blue")

lines(c(min.pos, max.pos), c(15, 15), lty = "dotted", lwd = 1, col = "blue")

axis(2, at = c(0, 10, 20), labels = c(0, 10, 20), las = 1, lwd = 2) 
# mtext("LOD", side = 2, at = 10, line = 2, font = 2)

close.screen(all.screens = TRUE)

dev.off()

# pdf("Pdf/GLM_SNP4.pdf", width = 9, height = 6)
#   make.locus.200k("marker", "4-119610606", gwas, hg19)
# dev.off()

# Gene: STAT4
# pdf("Pdf/STAT4_GLM.pdf", width = 9, height = 6)
#   make.locus.200k("gene", "STAT4", gwas, hg19)
# dev.off()
