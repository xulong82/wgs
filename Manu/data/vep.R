rm(list = ls())
setwd("~/Dropbox/GitHub/Adsp")
load("./Manu/data.rdt")

options(stringsAsFactors = F)

so <- read.delim("./data/so.txt")
so <- gsub(" ", "", so$SO.term)

vep <- data$vepAll

cons <- sapply(so, function(y) sapply(vep, function (x) sum(x[grep(y, names(x))])))
cons <- colSums(cons)
cons <- cons[! cons == 0]

sort(cons / sum(cons), decreasing = T)

pdf("./Pdf/cons.pdf", width = 12, height = 8)

op <- par(mar = c(5, 20, 4, 2))
bar <- barplot(cons, xlim = c(0, max(cons) + 5e6), 
               axes = F, border = NA, horiz = T, las = 1, space = 0.75)
abline(0, 1, lwd = 1, col = "black")
text(y = bar, x = cons + 1.8e6, labels = cons)

dev.off()

vep <- data$top$vep

cons <- sapply(so, function(y) sum(grepl(y, vep$Consequence)))
cons <- cons[! cons == 0]

sort(cons / sum(cons))

genes <- vep.lod$Symbol %>% unique
genes <- genes[genes != "-"]

pdf("./Manu/cons_glm.pdf", width = 8, height = 4.5)

op <- par(mar = c(5, 20, 4, 3))
bar <- barplot(cons, xlim = c(0, max(cons) + 1.5e2), 
               axes = F, border = NA, horiz = T, las = 1, space = 0.75)
abline(v = 0, lwd = 1, col = "black")
text(y = bar, x = cons + 5e1, labels = cons)

dev.off()

vep <- epis$epis2$vep

cons <- sapply(so, function(y) sum(grepl(y, vep$Consequence)))
cons <- cons[! cons == 0]

pdf("./Manu/cons_epis2.pdf", width = 8, height = 5)

op <- par(mar = c(5, 20, 4, 3))
bar <- barplot(cons, xlim = c(0, max(cons) + 1.5e2), 
               axes = F, border = NA, horiz = T, las = 1, space = 0.75)
abline(v = 0, lwd = 1, col = "black")
text(y = bar, x = cons + 7e1, labels = cons)

dev.off()

vep <- epis$epis4$vep

cons <- sapply(so, function(y) sum(grepl(y, vep$Consequence)))
cons <- cons[! cons == 0]

pdf("./Manu/cons_epis4.pdf", width = 8, height = 4)

op <- par(mar = c(5, 20, 4, 3))
bar <- barplot(cons, xlim = c(0, max(cons) + 1.5e2), 
               axes = F, border = NA, horiz = T, las = 1, space = 0.75)
abline(v = 0, lwd = 1, col = "black")
text(y = bar, x = cons + 1e1, labels = cons)

dev.off()

