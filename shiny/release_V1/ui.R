library(shiny)
library(ggvis)

load("so.rdt")

shinyUI(
  navbarPage("ADSPedia",
     tabPanel("Search a locus",
       fluidPage(
              p(strong("A"), "scatterplot of LOD in 5 Mbp genomic range. Variants with LOD values bigger than 6 are included. (red, SNP; green, IND; gene annotation, hg19 assembly)."),
              hr(),
              tags$div(textInput("gene", label = h5(""), value = "SORBS2"), align = "center", 
                       helpText("Please input a gene symbol or a genomic position (1-123456).")),
              plotOutput("loci")
     )),
     tabPanel("Explore the variants", 
       fluidPage(
         p(strong("An"), "interactive scatterplot of variants in efect size, MAF, and LOD. BiSNP and BiIND with LOD bigger than 10 are included for query."),
         hr(),
         fluidRow(
           column(3,
             checkboxInput("vep_all", "Select/deselect all"), 
             checkboxGroupInput("vep", "by effects", choices = so, selected = so)
           ),
           column(9,
              fluidRow(
                column(2, 
                  selectInput("x", "X-axis", choices = list("LOD" = "LOD", "Effect" = "PAR", "MAF" = "MAF"), selected = "MAF")
                ),
                column(2, 
                  selectInput("y", "Y-axis", choices = list("LOD" = "LOD", "Effect" = "PAR", "MAF" = "MAF"), selected = "PAR")
                ),
                column(4, 
                  sliderInput("lod", "Lod range", min = 10, max = 40, value = c(20, 40))
                ),
                column(4, 
                  sliderInput("par", "Effect size", min = 0, max = 4, value = 0.5, step = 0.1)
                )
              ),    
              hr(),
              ggvisOutput("ggvis"),
              br(),
              tableOutput("single_var")
           ))
         )),
     tabPanel("Explore the genes",
              fluidPage(
              p(strong("Relevant"), "variants of a gene. Please input a gene symbol to return its Entrez summary and relevant variants."),
              hr(),
              tags$div(textInput("gene1", label = h5(""), value = "SORBS2"), align = "center", helpText("All captalized, such as SORBS2")),
              hr(),
              fluidRow(
                column(4, tags$div(sliderInput("lod1", "Lod range", min = 10, max = 40, value = c(15, 40)), align = "center")),
                column(4, tags$div(sliderInput("par1", "Effect size", min = 0, max = 4, value = 0.2, step = 0.1), align = "center")),
                column(4, tags$div(align = "center", checkboxGroupInput("vep1", "Remove mutants by consequences?", 
                       choices = list("intron_variant" = "intron_variant", "intergenic_variant" = "intergenic_variant"))))
              ),
              hr(),
              textOutput("summary"),
              hr(),
              tableOutput("gene1")
     )),
     tabPanel("Manhattan",
              fluidPage(
                p(strong("Manhattan"), "plot of the GWAS with generalized mixed linear model. Variants with LOD values bigger than 10 are included."),
                helpText("Plan to implement interactive functionalities on this plot. There are technical issues on auto scaling with ggvis."),
                hr(),
                plotOutput("gwas")
              ))   
))
    
