library(ggvis)
library(shiny)

shinyUI(
  navbarPage("ADSPedia",
     tabPanel("Explore", 
         fluidPage(
           p(strong("Candidate"), "variants to make APOE4-sensitized mouse model for Alzheimer's disease."),
           hr(),
           fluidRow(
             column(2,
               selectInput("x", "X-axis", choices = list("LOD" = "LOD", "Effect" = "PAR", "MAF" = "MAF"), selected = "LOD"),
               selectInput("y", "Y-axis", choices = list("LOD" = "LOD", "Effect" = "PAR", "MAF" = "MAF"), selected = "PAR")
             ),
             column(10,
               ggvisOutput("ggvis"),
               br()
           )),
               tableOutput("single1"),
               tableOutput("single2")
     )),
     tabPanel("Gene",
              fluidPage(
                p(strong("Candidate"), "genes to make APOE4-sensitized mouse model for Alzheimer's disease."),
                hr(),
                tableOutput("summary1")
              )),
     tabPanel("Variant",
              fluidPage(
                p(strong("Candidate"), "variants to make APOE4-sensitized mouse model for Alzheimer's disease."),
                hr(),
                tableOutput("summary2")
              ))
))
    
