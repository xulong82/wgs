load("/data/xwang/Adsp/GWAS/optimizing.rdt")
load("/data/xwang/byglmm/igap.rdt")

gwas = do.call(rbind, glm.fit)
igap = igap_s1[match(gwas$UID, igap_s1$UID), ]

p0 = igap$Pvalue
p1 = pchisq(gwas$LOD, df = 1, lower.tail = F)

save(p0, p1, file = "/data/xwang/byglmm/weight.rdt") 

# adj = iGWAS(P_current = pval, N_current = 576, P_prior = igap$Pvalue, N_prior = N0)
