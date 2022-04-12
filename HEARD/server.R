library("shiny")
library("tidyverse")
library("readxl")
library("readr")
library("dplyr")
library("tools")
library("DT")
library("bslib")
library("shinycustomloader")

options(shiny.maxRequestSize=300*1024^2)
options(repos = BiocManager::repositories())

function(input, output, session) {

  uploaded_HRD_data <- eventReactive(c(input$uploaded_HRD_data),
                                         { if (file_ext(input$uploaded_HRD_data$datapath) == "xlsx") {
                                              data_viz <- read_excel(input$uploaded_HRD_data$datapath)
                                           
                                              data_viz
                                         }
                                           else {
                                             data_viz <- read_delim(input$uploaded_HRD_data$datapath, show_col_types = FALSE)
                                             
                                             data_viz
                                           }
                                         })
  
  
  output$dt_PatientData <- renderDataTable(server = TRUE,{
    input$anaylze_data
    uploaded_HRD_data <- isolate(uploaded_HRD_data())
    datatable(
      uploaded_HRD_data,
      extensions = c("Buttons", "Scroller"),
      filter = "top",
      fillContainer = TRUE,
      rownames = FALSE,
      options = list(
        searching = TRUE,
        deferRender = FALSE,
        scrollY = 200,
        scroller = TRUE,
        autoWidth = TRUE,
        ordering = TRUE,
        buttons = c("colvis", "copy", "csv", "excel", "print", "pdf"),
        dom = "Bfrtip"))
    })
}