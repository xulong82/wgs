library(dplyr)

rm(list = ls())
hpc <- "/data/xwang/Adsp"
github <- "~/Dropbox/GitHub/Adsp"

setwd(github)
load("data/glmList.rdt")
load("data/glmList_rare.rdt")

gwas_lod <- filter(glmList$gwas, LOD > 15) # permutation cut

var <- gwas_lod$UID

chr <- paste0("chr", 1:22)

setwd(hpc)

geno <- lapply(chr, function(chr) { 
  cat(chr, "\n")
  load(paste0("biVar/biGeno_", chr, ".rdt"))
  biGeno[rownames(biGeno) %in% var, ]
})

geno <- do.call(rbind, geno)

glmList$geno_top = geno

# --- variants by peaks ---

setwd(github)
load("data/glmList.rdt")

group <- glmList$group

peaks <- unique(as.character(group$PEAK))
peaks <- peaks[ ! grepl("X", peaks) ]

setwd(hpc)

geno <- sapply(peaks, function(x) {
  cat(x, "\n")

  chr = gsub("-.*", "", x)
  pos = gsub(".*-", "", x) %>% as.numeric
  limit1 = pos - 1e5
  limit2 = pos + 1e5

  load(paste0("golden/geno_chr", chr, ".rdt"))
  load(paste0("golden/meta_chr", chr, ".rdt"))

  UID = meta$UID[ (meta$POS > limit1) & (meta$POS < limit2) ]
  geno[UID, ]

})

save(geno, file = "./GWAS/geno.rdt")

