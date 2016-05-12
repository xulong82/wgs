library(shiny)
library(ggvis)
library(dplyr)
library(ggplot2)

load("data.rdt")
for(obj in names(shinyList)) assign(obj, shinyList[[obj]])

shinyServer(function(input, output, session) {
  
    table <- addi[c(17:22, 16, 7, 25:29, 33:35, 38, 39)]
    points <- reactiveValues(selected = NULL)
    
    output$table_1 <- renderTable({
      table_1 <- addi[c(17:22, 16, 39, 7, 25:29, 33:36, 38)]
      table_1$Consequence <- gsub(",", ", ", table_1$Consequence)
      table_1
    }, include.rownames = FALSE)
    
    output$table_2 <- renderTable({
      table_2 <- apoe4[c(17:22, 16, 39, 7, 25:29, 33:36, 38)]
      table_2$Consequence <- gsub(",", ", ", table_2$Consequence)
      table_2
    }, include.rownames = FALSE)
    
    output$table_3 <- renderTable({
      table_3 <- apoe2[c(17:22, 16, 39, 7, 25:29, 33:36, 38)]
      table_3$Consequence <- gsub(",", ", ", table_3$Consequence)
      table_3
    }, include.rownames = FALSE)
    
    tooltip_values <- function(x) {
      if (is.null(x)) return (NULL)
      points$selected <- x$UID
      entry <- table[table$UID == x$UID, c(1:4)]
      paste0(names(entry), ": ", format(entry), collapse = "<br />")
    }
    
    gv1 <- reactive({
      table %>%
        ggvis(~MAF, ~pSnp, key := ~UID) %>%
        layer_points() %>%
        add_tooltip(tooltip_values, "hover") 
    }) %>%  bind_shiny("ggvis")
    
    output$note1 <- renderTable({
      table[table$UID == points$selected, 1:17]
    }, include.rownames = FALSE)

    output$note2 <- renderTable({
      table[table$UID == points$selected, c(7, 18)]
    }, include.rownames = FALSE)

})
