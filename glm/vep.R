# VEP INPUT

load("/data/xwang/adsp3/sampling.rdt")
mcmc_id = gsub("_.*", "", rownames(mcmc)) # 9726

vep.input = read.table("/data/xwang/adsp3/plink/wgs.bim")
vep.input = vep.input[vep.input$V2 %in% mcmc_id, ]

vep.input$V7 = paste(vep.input$V5, vep.input$V6, sep = "/")
vep.input = vep.input[c("V1", "V4", "V4", "V7")]
vep.input$V8 = "+"

write.table(vep.input, file = "~/Dropbox/GitHub/wgs/vep/vep_input.txt", quote = F, sep = "\t", row.names = F)

vep = read.table("~/Dropbox/GitHub/wgs/vep/vep_output.txt", stringsAsFactors = F, header = T)

save(vep, file = "~/Dropbox/GitHub/wgs/vep/vep.rdt")

# GRAPH AND STATISTICS

