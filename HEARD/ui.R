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
options(shiny.maxRequestSize=300*1024^2)
options(repos = BiocManager::repositories())

tcga_data<-read.table(file="./TCGA_HRD_positive_samples.txt", sep="\t", header=TRUE)
sbs_score<-read.table(file="./TCGA_SBS_signature_exposure.txt", sep="\t", header=TRUE)
id_score<-read.table(file="./TCGA_ID_signature_exposures.txt", sep="\t", header=TRUE)
gene_loh <- read.table(file = "./gene_LOH_events.txt",header = TRUE)



#EVAN

clinical_data <- read.table(file="./TCGA_HRD_positive_samples_mut_calls.txt", sep="\t",header=TRUE)
hrd_data <- read.table(file="./TCGA_HRD_positive_samples.txt", sep="\t", header=TRUE)
id_data<-read.table(file="./TCGA_ID_signature_exposures.txt", sep="\t", header=TRUE)
sbs_data <- read.table(file="./TCGA_SBS_signature_exposure.txt", sep="\t", header=TRUE)
mut_call_data <- read.table(file="./TCGA_HRD_positive_samples_mut_calls.txt", sep="\t", header=TRUE)


#EVAN
# A couple of convenience utilities
chr_translation=list("-01","-02","-03","-04","-05",
                     "-06","-07","-08","-09","-10",
                     "-11","-12","-13","-14","-15",
                     "-16","-17","-18","-19","-20",
                     "-21","-22","-23")

names(chr_translation)=c("chr1","chr2","chr3","chr4","chr5",
                         "chr6","chr7","chr8","chr9","chr10",
                         "chr11","chr12","chr13","chr14","chr15",
                         "chr16","chr17","chr18","chr19","chr20",
                         "chr21","chr22","chrX")


#utility functions
getShortName<-function(sample_basename) { 
  tokens=str_split(sample_basename,"\\.",n=Inf)
  parta=unlist(tokens)[1]
  tokensb=str_split(parta,"-")
  tokensb=unlist(tokensb)
  shortname=paste(tokensb[1],tokensb[2],tokensb[3],sep="-")
  return(shortname)
}




navbarPage(
  theme = bs_theme(version = 5),
  "HEARD",
  img(
    src = "images/HEARD.png",
    width = "150px"
  ),
  tabPanel("Patient_View",
           sidebarPanel(width=2,
                        # naming this panel the Case Navigator
                        h3("Case Navigator"),
                        # get patient ID info for rendering
                        pickerInput("pat_id", "Patient:", 
                                    choices=(sbs_data$FileName)),
                        # Radio button for report
                        radioButtons("hrd_stat","HRD Positive",choices=c("Yes","No")),
                        # button for downloading the report
                        h3("Click to Generate Report"),
                        downloadButton("report", "Generate report"),
           ),
                                # First panel is the patient data - NEED TO RENDER
                                fluidRow(
                                  column(
                                    width=12,offset=0,
                                    #withLoader(textOutput("clinical_info")),
                                    htmlOutput("clinical_info")           
                                  ), # column end
                                ),
                                fluidRow(
                                  # First Row is the HRD scores
                                  column(width=3,
                                         h3("HRD Score Metrics"),
                                         pickerInput("hrd_metrics", "HRD Metric:", 
                                                     choices=colnames(hrd_data[,-1])),
                                         withLoader(plotOutput("hrdScores", height="400px", width="400px"),
                                                    type="html", loader="pacman")
                                  ), #end hrdScores column
                                  column(width=3,offset=3,
                                         h3("Mutational Signature Exposure"),
                                         withLoader(plotOutput("mutsigs", height="400px", width="400px"),
                                                    type="html", loader="pacman")
                                  ),
                                ),
                                # next Row is the chromosome info
                                fluidRow(
                                  column(width=2,offset=6,
                                         pickerInput("Chromosome", "Chromosome:",
                                                     choices = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13',
                                                                 'chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX'))),
                                ),
                                fluidRow(
                                  column(width=12, offset=0,
                                         withLoader(plotOutput("chromImage", height="300px", width="300px"),
                                                    type="html", loader="pacman")
                                  ), #end chromImage column
                                ),
                                #), #end fluidRow2
                                #),#end fluidrow1
                       ),
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

