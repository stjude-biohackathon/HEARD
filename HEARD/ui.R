library("shiny")
library("tidyverse")
library("readxl")
library("readr")
library("dplyr")
library("tools")
library("DT")
library("bslib")
library("shinycustomloader")
library("plotly")
library("shinyWidgets")
library("ggplot2")
library("RColorBrewer")
library("ComplexHeatmap")

options(repos = BiocManager::repositories())

tcga_data<-read.table(file="./TCGA_HRD_positive_samples.txt", sep="\t", header=TRUE)
sbs_score<-read.table(file="./TCGA_SBS_signature_exposure.txt", sep="\t", header=TRUE)
id_score<-read.table(file="./TCGA_ID_signature_exposures.txt", sep="\t", header=TRUE)
gene_loh <- read.table(file = "./gene_LOH_events.txt",header = TRUE)


navbarPage(
  theme = bs_theme(version = 5),
  "HEARD",
  img(
    src = "images/HEARD.png",
    width = "150px"
  ),
  tabPanel("Patient_View"),
  tabPanel("Cohort_View",
           sidebarPanel(pickerInput("TCGA", "TCGA:", 
                                    choices = (tcga_data[,1])),
                        pickerInput("Column", "Column:",
                                    choices=colnames(tcga_data)[!colnames(tcga_data) %in% c("FileName")])
                        
           ),
           #We can change the actual loader to be a gif or image
           fluidRow(
             column(
               width=12,
               withLoader(plotOutput("LOH_Info"), type="html", loader="pacman")
             )
           ),
           fluidRow(
             column(
               width=12,
               withLoader(plotOutput("scoreTCGA"), type="html", loader="dnaspin")
             )
           ),
           fluidRow(
             column(
               width=5,
               withLoader(plotOutput("id_htmap"), type="html", loader="dnaspin")
             ),
             column(
               width=7,
               withLoader(plotOutput("sbs_htmap"), type="html", loader="dnaspin")
             )
           )
  ),
  tabPanel(
    "Upload Data Table",
    fluidPage(
      wellPanel(p("Upload your data file below and click Analyze Data to begin visualization.")),
      fluidRow(
        column(fileInput("uploaded_HRD_data",
                         label = "Upload your data",
                         accept = c(".csv", ".tsv", ".txt", ".xlsx")),
               width = 6)
      ),
      actionButton("anaylze_data",
                   "Analyze Data"),
      withLoader(dataTableOutput("dt_PatientData", width = "100%", height = "auto"), type="html", loader="pacman")
    )
  ),
  
  collapsible = FALSE
  
  
)
