rm(list = ls())

setwd("~/Dropbox/GitHub/Adsp")
load("./data/mdata.rdt")
SRR <- mdata$SRR
Family <- mdata$Family.ID
Apoe2 <- mdata$Apoe2
Apoe3 <- mdata$Apoe3
Apoe4 <- mdata$Apoe4

chr <- commandArgs(TRUE)

setwd("/data/xwang/ADSP/biVar")
load(paste("biInf_chr", chr, ".rdt", sep = ""))
load(paste("biGeno_chr", chr, ".rdt", sep = ""))

biInf$SNP <- "N"
NC <- c("A", "T", "C", "G")
biInf$SNP[biInf$REF %in% NC & biInf$ALT %in% NC] <- "Y"

maf <- rowSums(biGeno) / ncol(biGeno) / 2
biInf$MA <- biInf$ALT
biInf$MA[maf > .5] <- biInf$REF[maf > .5]

maf[maf > .5] <- 1 - maf[maf > .5]
biInf$MAF <- maf
biInf$HET <- apply(biGeno, 1, function (x) sum(x == 1)) / ncol(biGeno)
  
biInf$N.SRR <- apply(biGeno, 1, function (x) length(SRR[x > 0]))
biInf$N.Family <- apply(biGeno, 1, function (x) length(unique(Family[x > 0])))

biInf$COR.APOE2 <- apply(biGeno, 1, function(x) cor(x, Apoe2))
biInf$COR.APOE3 <- apply(biGeno, 1, function(x) cor(x, Apoe3))
biInf$COR.APOE4 <- apply(biGeno, 1, function(x) cor(x, Apoe4))

meta <- biInf
save(meta, file = paste0("../META/meta_chr", chr, ".rdt"))

