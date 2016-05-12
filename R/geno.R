library(dplyr)

rm(list = ls())

setwd("/data/xwang")
load("./bgwas/igap.rdt")

var <- rownames(igap_s1) <- igap_s1$UID

chr <- paste0("chr", 1:22)

# shared variants of ADSP and IGAP, and prior
 
lapply(chr, function(chr) { cat(chr, "\n")
  load(paste0("./Adsp/golden/geno_", chr, ".rdt"))
  geno = geno[rownames(geno) %in% var, ]
  igap = igap_s1[match(rownames(geno), igap_s1$UID), c("Beta", "SE")]
  save(geno, file = paste0("./bgwas/geno/", chr, ".rdt"))
  save(igap, file = paste0("./bgwas/prior/", chr, ".rdt"))
})

