library(ggvis)
library(shiny)

shinyUI(
  navbarPage("ADSPedia",
     tabPanel("Explore", 
         fluidPage(
           p(strong("Candidate"), "variants to make mouse model for Alzheimer's disease."),
           hr(),
               ggvisOutput("ggvis"),
               br(),
               tableOutput("note1"),
               tableOutput("note2")
           )),
     tabPanel("Additive",
              fluidPage(
                p(strong("Top"), "variants (LOD value in 0.001 quantile) that affect protein coding in additive model."),
                hr(),
                tableOutput("table_1")
              )),
     tabPanel("epi_Apoe4",
              fluidPage(
                p(strong("Most"), "risky variants who also affect protein coding. LOD value in top 0.001 quantile."),
                hr(),
                tableOutput("table_2")
              )),
     tabPanel("epi_Apoe2",
              fluidPage(
                p(strong("Most"), "protective variants who also affect protein coding. LOD value in top 0.001 quantile."),
                hr(),
                tableOutput("table_3")
              ))
))
    
