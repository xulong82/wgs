library(shiny)
library(ggvis)
library(ggplot2)

load("hg19.rdt")
load("lod.rdt")
load("table.rdt")
load("summary.rdt")
load("so.rdt")

shinyServer(function(input, output, session) {
  
  table$cons <- rep("Protective", nrow(table))
  table$cons[sign(table$PAR) == 1] <- "Risky"
  
  points <- reactiveValues(selected = NULL)
  tooltip_values <- function(x) {
    if (is.null(x)) return (NULL)
    points$selected <- x$UID
    entry <- table[table$UID == x$UID, c(2:4, 9:10)]
    paste0(names(entry), ": ", format(entry), collapse = "<br />")
  }
  
  observe({
    updateCheckboxGroupInput(session, "vep", choices = so, selected = if (! input$vep_all) so)
  })
  
  graph.dt <- reactive({
    graph.dt <- table[table$LOD > input$lod[1] & table$LOD < input$lod[2], ]
    graph.dt <- graph.dt[abs(graph.dt$PAR) > input$par, ]
    graph.dt <- graph.dt[sapply(strsplit(graph.dt$Consequence, ","), function (x) any(x %in% input$vep)), ]
    graph.dt$x <- graph.dt[, input$x]
    graph.dt$y <- graph.dt[, input$y]
    graph.dt
  })
  
  gv1 <- reactive({
    graph.dt %>%
    ggvis(~x, ~y, shape = ~cons, fill = ~cons, key := ~UID) %>%
    layer_points() %>%
    add_axis("x", title = input$x, title_offset = 40) %>%
    add_axis("y", title = input$y, title_offset = 40) %>%
    add_legend(c("shape", "fill"), title = "") %>%
    add_tooltip(tooltip_values, "hover") %>%
    set_options(height = 400, width = 500) 
  }) %>%  bind_shiny("ggvis")
  
  glm <- var.df[var.df$LOD > 10, ]
  glm$POS <- glm$POS * 1e-3  # mb
  glm$CHR1 <- gsub("X", 23, glm$CHR)
  glm$CHR1 <- as.numeric(gsub("Y", 24, glm$CHR1))
  glm$POS1 <- c(0, chrlen)[glm$CHR1] + glm$POS
  glm$COL <- rep(0, nrow(glm))
  glm$COL[glm$CHR1 %% 2 == 1] <- 1
  mycol <- c("dodgerblue3", "firebrick1")
  
  output$gwas <- renderPlot({
    ggplot(glm, aes(x = POS1, y = LOD, color = as.factor(COL))) +
    geom_point(alpha = 0.7) + guides(color = F) +
    theme_bw() + xlab("") + ylab("LOD") +
    scale_color_manual(values = mycol) +
    scale_x_continuous(breaks = chrmid, labels = names(chrlen)) +
    theme(axis.text.x = element_text(color = rep(rev(mycol))))
  }, height = 500)
  
  output$single_var <- renderTable({
    if (is.null(points$selected)) 
      var_entry <- data.frame("Hove_over_a_point_for_details" = "")
    else
      var_entry <- table[table$UID == points$selected, -c(1, 16, 18)]
    var_entry
  }, include.rownames = FALSE)

  output$loci <- renderPlot({
    name <- input$gene
    if (grepl("-", name)) {
      chr <- as.numeric(gsub("-.*", "", name))
      pos <- as.numeric(gsub("^.*-", "", name)) * 1e-3
    } else {
      entry <- hg19[hg19$NAME == name, ]
      chr <- entry$CHR
      pos <- entry$POS
    }
    minp <- pos - 2.5e3
    maxp <- pos + 2.5e3
    
    var <- subset(var.df, var.df$CHR == chr & var.df$POS >= minp & var.df$POS <= maxp)
    var <- var[order(var$POS), ]
    gene <- subset(hg19, hg19$CHR == chr &
                     ((hg19$START >= minp & hg19$START <= maxp) | (hg19$END >= minp & hg19$END <= maxp)))
    gene <- gene[order(gene$POS), ]
    
    snp <- var[var$TYPE == "SNP", ]
    ind <- var[var$TYPE == "IND", ]
    
    plot(snp$POS, snp$LOD, type = "p", pch = 23, bg = "red", main = "", 
         xlab = "", ylab = "", xlim = c(minp, maxp), ylim = c(-30, 39), axes = F)
    box()
    points(ind$POS, ind$LOD, type = "p", pch = 23, bg = "green", main = "", 
         xlab = "", ylab = "", xlim = c(minp, maxp), ylim = c(-30, 39), axes = F)
    lines(c(minp, maxp), c(0, 0), lty = "dotted", lwd = 1, col = "black")
    lines(c(minp, maxp), c(24, 24), lty = "dotted", lwd = 1, col = "blue")
    axis(1, at = c(minp, pos, maxp), labels = round(c(minp, pos, maxp)), las = 1) 
    mtext(paste("Chromosome", chr, "Position (Kb)"), side = 1, line = 2.5, font = 2)
    axis(2, at = c(0, 6, 12, 18, 24, 30, 36), labels = c(0, 6, 12, 18, 24, 30, 36), las=1) 
    mtext("LOD", side = 2, at = 9, line = 2, font = 2)
    if (nrow(gene) != 0) {
      for (i in 1:nrow(gene)) {  # plot the genes
        adj.text <- -(i %% 15 + 1) * 1.7 
        if (!is.na(gene[i, ]$NAME))
          text(gene[i, ]$POS, adj.text, labels = gene[i, ]$NAME, cex = 0.8, pos = 1, offset = 0, font = 2)
      }
    }
  }, height = 600)
  
  output$summary <- renderText({
    paste(query[query$query == input$gene1, ], collapse = ", ")
  })

  output$gene1 <- renderTable({
    table1 <- table[table$Symbol == input$gene1, -c(1, 18)]
    table1 <- table1[! duplicated(table1), ]
    if ("intron_variant" %in% input$vep1) table1 <- table1[! grepl("intron_variant", table1$Consequence), ]
    if ("intergenic_variant" %in% input$vep1) table1 <- table1[! table1$Consequence == "intergenic_variant", ]
    table1 <- table1[table1$LOD > input$lod1[1] & table1$LOD < input$lod1[2], ]
    table1 <- table1[abs(table1$PAR) > input$par1, ]
    rownames(table1) <- NULL
    table1
  })

})

