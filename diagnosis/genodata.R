# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: generate genotype data for STAN

rm(list = ls())

#--- LogLik in single scan
load("/data/xwang/ADSP/R/sscan_clm_rs.rdt")  # rs-prefixed
sscan1 <- sscan
load("/data/xwang/ADSP/R/sscan_clm_novo.rdt")  # de-novo
sscan$CHROM <- paste("chr", sscan$CHROM, sep = "")
sscan <- rbind(sscan1, sscan[, c("CHROM", "POS", "ID", "m1", "m2", "m3", "m4")])
rm(sscan1)

snpId <- "rs34572242"

entry <- sscan[sscan$ID == snpId, ]
chromosome <- entry$CHROM 
position <- entry$POS
snpIdPool <- subset(sscan,  # +/- 10m
  sscan$CHROM == chromosome & sscan$POS >= position - 1e5 & sscan$POS <= position + 1e5)$ID

load("/data/xwang/ADSP/R/vdata_rs.rdt")
geno1 <- geno[rownames(geno) %in% snpIdPool, ]
load("/data/xwang/ADSP/R/vdata_denovo.rdt")
geno <- rbind(geno1, geno[rownames(geno) %in% snpIdPool, ])
rm(geno1)

save(geno, file = "~/Dropbox/ADSP/Stan/genoSNPrs34572242.rdt")

