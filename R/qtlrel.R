library(QTLRel)

setwd("~/GitHub/Adsp/LMM")

rm(list = ls())
load("miscEx.RData")

id = rownames(pdatF8)

idcf = cic(pedF8, ids = id, df = 3, ask = T, verbose = T)

names(gmF8)

idx = ! is.na(pdatF8[, "bwt"])
pdatTmp = pdatF8[idx, ]

gdatTmp<- gdatF8[match(rownames(pdatTmp),rownames(gdatF8)), ]

ii<- match(rownames(pdatTmp),rownames(gmF8$AA))

vc<- estVC(y = pdatTmp[,"bwt"], x = pdatTmp[,c("sex","age")],
           v=list(AA=gmF8$AA[ii,ii],DD=gmF8$DD[ii,ii],HH=NULL,AD=NULL,MH=NULL,
           EE=diag(nrow(pdatTmp))))

sum(is.na(gdatTmp))

gdatTmpImputed<- genoImpute(gdatTmp,gmapF8,gr=8,na.str=NA)

lrt<- scanOne(y=pdatTmp[,"bwt"], x=pdatTmp[,c("sex","age")],
              gdat=gdatTmpImputed, vc=vc)

plot(lrt,gmap=gmapF8,main="Body Weight") # plotting

plot(lrt$p,gmap=gmapF8,main="Body Weight") # plotting

lrt1 = lrt
plot(lrt1,gmap=gmapF8,main="Body Weight") # plotting

library(DOQTL)
library(MUGAExampleData)

rm(list = ls())
data(pheno)
data(model.probs)

K = kinship.probs(model.probs)

covar = data.frame(sex = as.numeric(pheno$Sex == "M"), diet = as.numeric(pheno$Diet == "hf"))
rownames(covar) = rownames(pheno)

load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))

qtl = scanone(pheno = pheno, pheno.col = "HDW2", probs = model.probs, K = K,
              addcovar = covar, snps = muga_snps)

plot(qtl, main = "HDW2")

lod = qtl$lod
