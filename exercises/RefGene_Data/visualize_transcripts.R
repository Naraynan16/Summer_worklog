
introns <- read.csv(file = "refgene_introns.csv",header = T)
utr <- read.csv(file = "all_utr.csv",header = T)
exons <- read.csv(file = "individual_transcripts.csv",header = T)

library(shiny)
library(shinythemes)
library(dplyr)
library(tidyr)

ui <- fluidPage(title = "Genomic Locations of Genes",
                sidebarLayout(
                  sidebarPanel(fluidRow(
                    column(5,
                           
                           textInput(inputId = "gene",
                                     label = "Enter gene name",
                                     value = "A1BG"),
                           
                           actionButton(inputId = "submit",
                                        label = "submit")
                           
              ))),
  
          mainPanel(
              tabsetPanel(type = "tabs",
                        tabPanel(dataTableOutput(outputId = "Exons"),title ="Exons"),
                        tabPanel(dataTableOutput(outputId = "Introns"),title = "Introns"),
                        tabPanel(dataTableOutput(outputId = "UTR"),title = "UTR")  
                        
                        ))))
      
          
server <- function(input, output)
  {
  
  Ex  <- eventReactive(input$submit,
                {
                  exons %>% filter(exons$geneName == input$gene)
                 } )
  In  <- eventReactive(input$submit,
                       {
                         introns %>% filter(introns$geneName == input$gene) 
                       } )
  Ut <- eventReactive(input$submit,
                      {
                         utr %>% filter(utr$geneName == input$gene)
                         
                         } )
  
  output$Exons <- renderDataTable(Ex())
  output$Introns <- renderDataTable(In())
  output$UTR <- renderDataTable((Ut()))
}
shinyApp(ui = ui, server = server)