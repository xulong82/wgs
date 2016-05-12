library(ggplot2)
library(GenomicRanges)
library(biomaRt)

rm(list = ls())
setwd("~/Dropbox/GitHub/Adsp")
load("Manu/data.rdt")

setwd("/data/xwang/Adsp")

chr <- c(1:22, "X", "Y")

gMix <- sapply(chr, function(x) { cat(x, "\n")
  load(paste0("gMix/gMix_chr", x, ".rdt"))
  sapply(gMix, nrow)
})

rownames(gMix) <- c("bi-SNP", "poly-SNP", "bi-Indel", "poly-Indel")

# ---

setwd("~/Dropbox/GitHub/Adsp")
gMix = data$gMix

gdt <- data.frame(num = c(gMix), type = rep(type, 24), chr = rep(chr, each = 4))

gdt$num <- gdt$num * 1e-6
gdt$chr <- factor(gdt$chr, levels = chr)
gdt$type <- factor(gdt$type, levels = type)

mycol <- c("grey70", "firebrick1", "dodgerblue3", "chartreuse3")
pdf("./Pdf/gMix_summary.pdf", width = 8, height = 4, family = "Helvetica")
ggplot(gdt, aes(x = chr, y = num, fill = type)) + 
  geom_bar(stat = "identity", width = 0.65) +
  theme_bw() + xlab("") + ylab("Variant Number in Million") +
  scale_fill_manual(values = mycol) +
  theme(panel.border = element_rect(size = 1, color = "grey30"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, vjust = 1),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

setwd("/data/xwang/Adsp")
mafInf <- sapply(chr, function(x) { cat(x, "\n")
  load(paste0("Meta/meta_chr", x, ".rdt"))
  snp = table(meta$MAF[meta$SNP == "Y"] < 0.01)
  ind = table(meta$MAF[meta$SNP == "N"] < 0.01)
  c(snp, ind)
})

rownames(mafInf) <- c("SNP > 0.01", "SNP < 0.01", "INDEL > 0.01", "INDEL < 0.01")

mafInf = data$mafInf
type = rownames(mafInf)

gdt <- data.frame(num = c(mafInf), type = rep(type, 24), chr = rep(chr, each = 4))

gdt$num <- gdt$num * 1e-6
gdt$chr <- factor(gdt$chr, levels = chr)
gdt$type2 <- factor(gsub("(SNP|INDEL).*", "\\1", gdt$type), levels = c("SNP", "INDEL"))
gdt$type3 <- factor(gsub(".*(> 0.01|< 0.01).*", "\\1", gdt$type), levels = c("< 0.01", "> 0.01"))

gdt1 <- gdt[gdt$type2 == "SNP", ]
gdt2 <- gdt[gdt$type2 == "INDEL", ]

setwd("~/Dropbox/GitHub/Adsp")

pdf("./Manu/snp_summary.pdf", width = 8, height = 4, family = "Helvetica")
ggplot(gdt1, aes(x = chr, y = num, fill = type3)) + 

pdf("./Manu/ind_summary.pdf", width = 8, height = 4, family = "Helvetica")
ggplot(gdt2, aes(x = chr, y = num, fill = type3)) + 

  geom_bar(stat = "identity", width = 0.65) + 
  theme_bw() + xlab("Chromosome") + ylab("Variant Number in Million") +
  scale_fill_manual(values = c("dodgerblue3", "firebrick1")) +
  theme(panel.border = element_rect(size = 1, color = "grey30"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), legend.key = element_blank()) 

dev.off()

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attr_gene = c("external_gene_name", "chromosome_name", "start_position", "end_position", "strand")
gene_ensembl = getBM(attr_gene, filters = "chromosome_name", values = c(1:22, "X", "Y"), mart = human)
attr_exon = c(att_gene, "exon_chrom_start", "exon_chrom_end")
exon_ensembl = getBM(attr_exon, filters = "chromosome_name", values = c(1:22, "X", "Y"), mart = human)

gene_gr = with(exon_ensembl, GRanges(chromosome_name, IRanges(start=start_position, end=end_position)))
(gene_gr_total = sum(width(disjoin(gene_gr))))
exon_gr = with(exon_ensembl, GRanges(chromosome_name, IRanges(start=exon_chrom_start, end=exon_chrom_end)))
(exon_gr_total = sum(width(disjoin(exon_gr))))
(intron_gr_total = gene_gr_total - exon_gr_total)
(genome_length = width(range(gene_gr)) %>% as.numeric %>% sum)
(intergenic_gr_total = genome_length - gene_gr_total)

exon_var = unlist(data$vep_UID_exon)
intron_var = unlist(data$vep_UID_intron)
intron_var = setdiff(intron_var, exon_var)

(exon_var_length = length(exon_var))
(intron_var_length = length(intron_var))
(intergenic_var_length = 12.6e6 - exon_var_length - intron_var_length)

(intergenic_density = intergenic_var_length / intergenic_gr_total)
(intron_density = intron_var_length / intron_gr_total)
(exon_density = exon_var_length / exon_gr_total)

