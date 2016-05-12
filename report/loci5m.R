make.locus <- function(type, name, marker.df, gene.df) {
  #--- extract locus
  if (type == "gene") {
    entry <- gene.df[gene.df$NAME == name, ]
  } else if (type == "marker") {
    entry <- marker.df[marker.df$ID == name, ]
  } else {
    stop("type = marker | name")
  }
  chromosome <- entry$CHR
  position <- entry$POS

  #--- locus in kb unit ---
  position <- position * 1e-3
  marker.df$POS <- marker.df$POS * 1e-3
  gene.df$START <- gene.df$START * 1e-3
  gene.df$END <- gene.df$END * 1e-3
  gene.df$POS <- gene.df$POS * 1e-3
  
  min.pos <- position - 2.5e3
  max.pos <- position + 2.5e3 

  #--- marker and genes in the locus ---
  marker <- subset(marker.df, marker.df$CHR == chromosome &
                   marker.df$POS >= min.pos & 
		   marker.df$POS <= max.pos)
  marker <- marker[order(marker$POS), ]
  
  gene <- subset(gene.df, gene.df$CHR == chromosome &
                          ((gene.df$START >= min.pos & gene.df$START <= max.pos) | 
                          (gene.df$END >= min.pos & gene.df$END <= max.pos)))
  gene <- gene[order(gene$POS), ]
  
  hit <- marker[marker$LOD == max(marker$LOD), ]
  if (nrow(hit) > 1) hit <- hit[1, ]
  
  #--- graphing ---
  par(mar = c(4, 4, 3, 4))
  # marker points
  plot(marker$POS, marker$LOD, type = "p", pch = 23, cex = 1.0, bg = "red", 
       main = "", xlab = "", ylab = "", xlim = c(min.pos, max.pos), ylim = c(-30, 39), axes = F)
  # box and lines
  box()
  lines(c(min.pos, max.pos), c(0, 0), lty = "dotted", lwd = 1, col = "black")
  lines(c(min.pos, max.pos), c(25, 25), lty = "dotted", lwd = 1, col = "blue")
  # top marker
  if (hit$LOD > 25) {
    points(hit$POS, hit$LOD, pch = 5, cex = 2.0, lwd = 2.0, col = "blue")
    text(hit$POS, hit$LOD, labels = hit$ID, pos = 3, offset = 1, font = 2)
  }
  # axis
  axis(1, at = c(min.pos, position, max.pos), 
       labels = round(c(min.pos, position, max.pos)), las = 1) 
  mtext(paste("Chromosome", chromosome, "Position (Kb)", sep=" "), side = 1, line = 2.5, font = 2)
  axis(2, at = c(0, 6, 12, 18, 24, 30, 36), labels = c(0, 6, 12, 18, 24, 30, 36), las=1) 
  mtext("LOD", side = 2, at = 9, line = 2, font = 2)
  # genes
  if (nrow(gene) != 0) {
    for (i in 1:nrow(gene)) {  # plot the genes
      adj.text <- -(i %% 15 + 1) * 1.7 
      if (!is.na(gene[i, ]$NAME))
        text(gene[i, ]$POS, adj.text, labels = gene[i, ]$NAME, cex = .5, pos = 1, offset = 0, font = 2)
    }
  }
}
