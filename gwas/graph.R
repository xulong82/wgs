rm(list = ls())

chrlen <- read.delim("~/Dropbox/GitHub/X/genomes/human.hg19.genome", header = F)
chrlen <- chrlen[match(paste0("chr", 1:22), chrlen$V1), ]
chrlen <- cumsum(as.numeric(chrlen$V2)) * 1e-6; names(chrlen) <- c(1:22)
chrmid <- diff(c(0, chrlen)) * 0.5 + c(0, chrlen[-length(chrlen)])

load("/data/xwang/Adsp/snptest.rdt")
load("/data/xwang/Adsp/GWAS/optimizing.rdt")

pval = pchisq(glm.fit, df = 1, lower.tail = F)

vId <- snptest$rsid
chr <- as.numeric(gsub("-.*", "", vId))
pos <- as.numeric(gsub(".*-", "", vId)) * 1e-6 # Mb
pos <- c(0, chrlen)[chr] + pos

pval = snptest$frequentist_add_pvalue
manhattan <- data.frame(uid = vId, chr = chr, pos = pos, pval = -log10(pval))
manhattan$col <- rep("a", nrow(manhattan)) 
manhattan$col[chr %% 2 == 1] <- "b"

png("~/Dropbox/GitHub/glmm/gwas/manhattan_snptest.png", width = 2e3, height = 1e3, res = 200)

ggplot(manhattan, aes(x = pos, y = pval, color = col)) + 
  geom_point(alpha = 0.9) + # geom_hline(yintercept = 7.3, color = "black") +
  scale_x_continuous(breaks = chrmid, labels = names(chrlen), limits = c(0, max(pos))) + 
  scale_color_manual(values = c("dodgerblue3", "firebrick1")) + 
  theme_bw() + xlab("") + ylab("-log10(P)") + guides(color = F) +
  theme(legend.key = element_blank())

dev.off()

