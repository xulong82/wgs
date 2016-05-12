library(dplyr)
library(ggplot2)
library(Biobase)

rm(list = ls())
hpc <- "/data/xwang/Adsp"
github <- "~/Dropbox/GitHub/Adsp"

setwd(hpc)

load("./GWAS/optimizing.rdt")
gwas <- do.call(rbind, glm.fit)

gwas <- gwas[, c("CHR", "POS", "SNP", "beta[5]", "LOD")]

setwd(github)

chr <- read.delim("../../X/genomes/human.hg19.genome", header = F)
chr <- chr[match(paste0("chr", 1:22), chr$V1), ]
# chr <- chr[match(paste0("chr", c(1:22, "X", "Y")), chr$V1), ]
chrlen <- cumsum(as.numeric(chr$V2)) * 1e-6  # Mb
names(chrlen) <- c(1:22)
# names(chrlen) <- c(1:22, "X", "Y")
chrmid <- diff(c(0, chrlen)) * 0.5 + c(0, chrlen[-length(chrlen)])

gwas$CHR1 <- gsub("X", 23, gwas$CHR) 
gwas$CHR1 <- as.numeric(gsub("Y", 24, gwas$CHR1))
gwas$POS <- gwas$POS * 1e-6 
gwas$POS1 <- c(0, chrlen)[gwas$CHR1] + gwas$POS

gwas$COL <- rep("Protective", nrow(gwas)) 
gwas$COL[gwas$"beta[5]" > 0] <- "Risky"

gwas$Effect = abs(gwas$"beta[5]")

gwas$Type = "SNP"
gwas$Type[gwas$SNP == "N"] = "INDEL"
gwas$Type = factor(gwas$Type, levels = c("SNP", "INDEL"))

mycol <- c("dodgerblue3", "firebrick1")

setwd(github)

date()
png("./Manu/manhattan.png", width = 2000, height = 1000, res = 200)

ggplot(gwas, aes(x = POS1, y = LOD, color = COL, shape = Type)) + 
  geom_point(alpha = 0.7) + 
# geom_point(alpha = 0.7, aes(size = Effect)) + 
  scale_x_continuous(breaks = chrmid, labels = names(chrlen)) +
  scale_color_manual(values = mycol) +
# guides(color = F, size = F) +
  theme_bw() + xlab("") + ylab("LOD") +
  theme(legend.title = element_blank(), 
        legend.key = element_blank())
dev.off()
date()

setwd(github)
load("data/glmList.rdt"); list <- glmList
for(obj in names(list)) assign(obj, list[[obj]])

# MANHATTAN PLOT

gwas <- filter(glm, LOD > quantile(LOD, 0.999))
gwas$pSnp <- lmm$snp  # QTLRel

gwas <- filter(hlm, LMER_LOD > quantile(LMER_LOD, 0.999))
gwas$LOD <- gwas$LMER_LOD
gwas$pSnp <- gwas$LMER_SIZE

gwas <- filter(hlm, CLMM_LOD > quantile(CLMM_LOD, 0.999))
gwas$LOD <- gwas$CLMM_LOD
gwas$pSnp <- gwas$CLMM_SIZE

gwas.lod <- filter(gwas, LOD > 15)
vep.lod <- filter(vep, UID %in% gwas.lod$UID)

genes <- vep.lod$Symbol %>% unique
genes <- genes[! genes == "-"]

variants <- sapply(genes, function(x) {
  vep_select = vep[vep$Symbol == x, ]
  vep_select$LOD = gwas$LOD[match(vep_select$UID, gwas$UID)]
  vep_select$UID[which.max(vep_select$LOD)]
})

variants["RP1-130L23"] <- ""
variants["CAMK2G"] <- ""
names(variants) <- gsub("NANOS1", "NANOS1\\*", names(variants))
gwas$LAB <- names(variants)[match(gwas$UID, variants)]
gwas$LAB[is.na(gwas$LAB)] <- ""

setwd(hpc)  
load("./GWAS/GLM.rdt") 
setwd(github) # GWAS peaks
source("./GWAS/peaks.R")
cutoff <- 15
chr <- paste0("chr", c(1:22, "X", "Y"))
peakList <- lapply(chr, function(x) { cat(x, "\n") # CWAS
  fit0 = glm.fit[[x]]
  hit0 = peak.detection(fit0$POS, fit0$LOD, cutoff)
  file.png <- paste0("./png/", x, ".png")
  png(file.png, width = 3000, height = 1500, res = 200)
  plot(fit0$POS, fit0$LOD, main = "", xlab = x, ylab = "LOD",
       ylim = c(0, 28), type = "p", pch = 19, cex = .5, col = "lightblue")
  if (length(hit0) > 0) {
    points(fit0$POS[hit0], fit0$LOD[hit0], type = "p", pch = 2, cex = 2, col = "red")
    text(fit0$POS[hit0], fit0$LOD[hit0], labels = fit0$UID[hit0], pos = 3, offset = 1)
  }; dev.off() 
  peak0 = fit0[hit0, ]
  min.pos = peak0$POS - 25e4  # 250K
  max.pos = peak0$POS + 25e4  # 250K
  cutoff = 3/4 * peak0$LOD
  ext0 <- do.call(rbind, lapply(1:nrow(peak0), function(x) {
    filter(fit0, POS > min.pos[x] & POS < max.pos[x] & LOD > cutoff[x])
  })); list(peak = peak0, ext = ext0)
}) 

