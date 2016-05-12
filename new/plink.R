library(matrixcalc)

setwd("/data/xwang/adsp3/plink")

date()
kin = read.table("wgs_kin.mibs")
date()

id = read.table("wgs_kin.mibs.id")

rownames(kin) = colnames(kin) = id$V1

is.positive.definite(as.matrix(kin))

save(kin, file = "/data/xwang/adsp3/plink/kin.rdt")

