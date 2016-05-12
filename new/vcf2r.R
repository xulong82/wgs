library(VariantAnnotation)

setwd("/data/xwang/adsp3/golden")

for ( chr in 1:8 ) {
  cat(paste0("chr", chr), "\n")
  fl = paste0("chr", chr, ".recode.vcf")
  geno = readGT(fl)
  geno[ geno == "0/0" | geno == "." ] = 0
  geno[ geno == "0/1" ] = 1
  geno[ geno == "1/1" ] = 2
  class(geno) = "numeric"
  save(geno, file = paste0("../R/chr", chr, ".rdt"))
}

