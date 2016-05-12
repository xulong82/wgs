library(shiny)
library(ggvis)
library(dplyr)
library(ggplot2)

load("beta.rdt")

shinyServer(function(input, output, session) {
  
    output$summary1 <- renderTable({
      beta$summary
    }, include.rownames = FALSE)
    
    output$summary2 <- renderTable({
      beta$variant[, -c(1, 10, 16, 24)]
    }, include.rownames = FALSE)
    
    table <- beta$variant
    points <- reactiveValues(selected = NULL)
    
    tooltip_values <- function(x) {
      if (is.null(x)) return (NULL)
      points$selected <- x$UID
      entry <- table[table$UID == x$UID, c(2:4)]
      paste0(names(entry), ": ", format(entry), collapse = "<br />")
    }
    
    ggvis.dt <- reactive({
      vis <- data.frame(UID = table$UID, LOD = table$LOD_epis, MAF = table$MAF,
        PAR = table$pSnp_epis + table$pApoe4_epis + table$pEpis_epis)
      vis <- vis[! duplicated(vis$UID), ]
      vis$x <- vis[, input$x]; vis$y <- vis[, input$y]
      vis
    })
    
    gv1 <- reactive({
      ggvis.dt %>%
        ggvis(~x, ~y, key := ~UID) %>%
        layer_points() %>%
        add_axis("x", title = input$x) %>%
        add_axis("y", title = input$y) %>%
        add_tooltip(tooltip_values, "hover") 
    }) %>%  bind_shiny("ggvis")
    
    output$single1 <- renderTable({
      if (is.null(points$selected)) 
        entry <- NULL
      else
        entry <- table[table$UID == points$selected, -c(1, 10, 16, 24)]
      entry
    }, include.rownames = FALSE)

    output$single2 <- renderTable({
      if (is.null(points$selected)) 
        entry <- NULL
      else {
        symbol <- table[table$UID == points$selected, "Symbol"]
        entry <- beta$summary[beta$summary$query == symbol, ]
      }
      entry
    }, include.rownames = FALSE)

})
