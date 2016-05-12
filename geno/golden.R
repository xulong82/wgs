rm(list = ls())

setwd("/data/xwang/Adsp")

chr <- paste0("chr", 1:22)

for (chr in chr) {

  cat(chr, "\n")

  load(paste0("Meta/meta_", chr, ".rdt"))
  load(paste0("biVar/biGeno_", chr, ".rdt"))
  
  all = all(rownames(biGeno) == meta$UID)
  cat(all, "\n")
  
  index <- which(meta$MAF > 0.01 & meta$HET < 0.99)
  
  meta = meta[index, ]
  geno = biGeno[index, ]
  
  save(meta, file = paste0("golden/meta_", chr, ".rdt"))
  save(geno, file = paste0("golden/geno_", chr, ".rdt"))

}

