library("shiny")
library("tidyverse")
library("readxl")
library("readr")
library("dplyr")
library("tools")
library("DT")
library("bslib")
library("shinycustomloader")

navbarPage(
  theme = bs_theme(version = 5),
  "HEARD",
  tabPanel(
    "Upload Data Table",
    fluidPage(
      img(
        src = "images/HEARD.png",
        width = "150px"
      ),
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