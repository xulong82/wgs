load("~/Dropbox/GitHub/wgs/loci/hg19.rdt") # HG19

make.locus.200k <- function(name, data) {

  # Locus in Kb
  data$P = -log10(data$P)
  data$POS = data$POS * 1e-3
  chr = data[data$UID == name, "CHR"]
  pos = data[data$UID == name, "POS"]
  min.pos = pos - 100; max.pos = pos + 100 
  
  # Markers and Genes
  md = subset(data, data$CHR == chr & data$POS >= min.pos & data$POS <= max.pos)
  gd = subset(hg19, hg19$CHR == chr & ((hg19$START >= min.pos & hg19$START <= max.pos) | (hg19$END >= min.pos & hg19$END <= max.pos)))
  gd = gd[! (gd$END > max.pos | gd$START < min.pos), ]
  gd = gd[order(gd$POS), ]
  hit = md[md$UID == name, ]
  
  # Plot
    par(mar = c(4, 4, 3, 4))
    plot(md$POS, md$P, type = "p", pch = 23, cex = 1.5, bg = "red", 
         main = "", xlab = "", ylab = "", xlim = c(min.pos, max.pos), ylim = c(-10, 10), axes = F)
    box()
    lines(c(min.pos, max.pos), c(0, 0), lty = "dotted", lwd = 1, col = "black")
    points(pos, md[md$UID == name, "P"], pch = 5, cex = 2.5, lwd = 2.5, col = "blue")
    text(pos, md[md$UID == name, "P"], label = name, pos = 3, offset = 1, font = 2)
    axis(1, at = c(min.pos, pos, max.pos), labels = round(c(min.pos, pos, max.pos)), las = 1) 
    mtext(paste0("Chromosome ", chr, " (Kb)"), side = 1, line = 2.5, font = 2)
    axis(2, at = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10), las=1) 
    mtext("-log10(P)", side = 2, at = 5, line = 2, font = 2)

    if (nrow(gd) != 0) {
      for (i in 1:nrow(gd)) {
        adj.arrow = -(i %% 3 + 1) * 2 
        adj.text = -(i %% 3 + 1) * 2 - 0.5 
        if (gd[i, ]$END - gd[i, ]$START > 1) {
          if (gd[i, ]$STRAND == "+" ) {
            arrows(max(gd[i, ]$START, min.pos), adj.arrow, min(gd[i, ]$END, max.pos), adj.arrow, 
                   length = 0.05, lwd = 2, code = 2, lty = "solid", col = "darkgreen")
          } else {    
            arrows(max(gd[i, ]$START, min.pos), adj.arrow, min(gd[i, ]$END, max.pos), adj.arrow, 
                   length = 0.05, lwd = 2, code = 1, lty = "solid", col = "darkgreen")
          }
        } else {
          points(gd[i, ]$POS, adj.arrow, pch = 16, col = "darkgreen")
        }
        
        if (!is.na(gd[i, ]$NAME)) {
          text(gd[i, ]$POS, adj.text, labels = gd[i, ]$NAME, cex = .6, pos = 1, offset = 0, font = 2)
        }
      }
    }  
}

