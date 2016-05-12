library(QTLRel)
library(parallel)

rm(list = ls())
hpc <- "/data/xwang/Adsp"
github <- "~/Dropbox/GitHub/Adsp"

arg <- commandArgs(TRUE) 
chr <- paste("chr", arg, sep = "")

setwd(github)
load("./data/mdata.rdt")
load("./data/kinship.rdt") 

KS <- kinship$autosome
if (arg %in% 1:22) KS <- kinship[[paste0("no_chr", arg)]]
KS[KS < 0] <- 0

setwd(hpc)
load(paste("biVar/biGeno_", chr, ".rdt", sep = "")) 
load(paste("Meta/meta_", chr, ".rdt", sep = "")) 

# uid <- meta$UID[meta$MAF > 0.01 & meta$HET < 0.99]
uid <- meta$UID[meta$MAF <= 0.01 | meta$HET >= 0.99]
geno <- biGeno[uid, ] + 1 # QTLRel use 1, 2, 3

n.core <- 15 
u.core <- round(nrow(geno) / n.core)
idx1 <- ((1:n.core) - 1) * u.core + 1 
idx2 <- c((1:(n.core - 1)) * u.core, nrow(geno))
gList <- mclapply(1:15, function(x) geno[idx1[x]:idx2[x], ], mc.preschedule = F)

# VC
Y <- mdata$AD1
X <- mdata[c("Age", "Sex", "Apoe2", "Apoe4")]
o <- estVC(y=Y, x=X, v=list(AA=KS, DD=NULL, HH=NULL, AD=NULL, MH=NULL, EE=diag(576)))

y <- mclapply(gList, mc.cores = n.core, FUN = function (gdat) {
       scanOne(Y, x = X, gdat = t(gdat), numGeno = TRUE, vc = o)})

assign(chr, y)
# save(list = chr, file = paste0("QTLRel/", chr, ".rdt"))
save(list = chr, file = paste0("QTLRel_rare/", chr, ".rdt"))

