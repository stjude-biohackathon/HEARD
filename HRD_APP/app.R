#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Goals- 
# Bar graph with interactive sorting (upgrade color selection)
# Gene List Frequency
# PNG incorporation-YES
# Interactive PCA plot with color selection of specific points
# Table with export and filtering<-Export and column filter, no search or sort
# Second tab with more information<-YES
#

library(shiny)
library(plotly)
library(shinyWidgets)
library(DT)
library(ggplot2)
library(shinycustomloader) #can include custom things during loading of elements
library(iBET)
options(repos = BiocManager::repositories())

patient_data<-read.table(file="./Sample_Excel_Shiny.txt", header=TRUE)
mutation_count<-read.table(file="./MutationCounts.txt", sep="\t")
tcga_data<-read.table(file="./TCGA_HRD_positive_samples.txt", sep="\t", header=TRUE)



#Create a function to make the shiny app

hrd_shiny<-function(patient_data) {
    body <- mainPanel(
        width = 10,
        fluidRow(
            column(width = 8,
                   withLoader(plotOutput("scorePlot", height = "500px", width = "500px"),
                              type = "html", loader = "dnaspin")
            ),
            column(width=8, 
                   withLoader(plotOutput("mutationCount", height="500px", width="500px"),
                               type="html", loader="pacman")
                
            ),
            column(width = 8,
                   #withLoader(plotOutput("dim_meta", height = "500px", width = "500px"),
                             # type = "html", loader = "dnaspin")
            ),
            column(width = 4,
                  # withLoader(plotOutput("dplot", height = "500px", width = "500px"),
                              #type = "html", loader = "dnaspin")
            )
        )
    )

# Define UI for application that draws a histogram


 ui <- fluidPage(

    # Application title
    titlePanel("Patient_Scores"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(width=2,
            h3("Plot Controls"),
            hr(style="margin:2px; background-color: #737373;"),
            pickerInput("Patient", "Patient:", 
                        choices = (patient_data[,1])),
            pickerInput("TCGA", "TCGA:", 
                      choices = (tcga_data[,1]))
    
        ),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(type="tabs",
                        tabPanel("Patient_Data",
                                 fluidRow(
                                     column(
                                         width=12,
                                         #div(h3("Patient_Data"), align="center"),
                                         #div(DT::dataTableOutput("dt_PatientData"), style="font-size: 100%; width: 100%")
                                         
                                     )
                                 )
                        ),
                        tabPanel("Cohort_View",
            #We can change the actual loader to be a gif or image
                            fluidRow(
                                column(
                                    width=12,
                                    withLoader(plotOutput("mutationCount"), type="html", loader="pacman")
                                    
                                )
                            ), 
                            fluidRow(
                                column(
                                    width=12,
                                    withLoader(plotOutput("scorePlot"), type="html", loader="dnaspin")
                                        ),
                                column(
                                    width=4,
                                    img(src='Tumor_Image.png', height=140, width=400)
                                        ),       
                                    )
                                ),
                fluidRow(
                    column(
                        width=12,
                        withLoader(plotOutput("scoreTCGA"), type="html", loader="dnaspin")
                    ),
                ),
                            tabPanel("Data_Tables",
                                     fluidRow(
                                         column(
                                             width=12,
                                             div(h3("Patient_Data"), align="center"),
                                             div(DT::dataTableOutput("dt_PatientData"), style="font-size: 100%; width: 100%")
                                             
                                             )
                                        )
                                     )
                                 )
                        )
                    )
)

# Define server logic required to draw a histogram
    server <- function(input, output, session) {

       
        
        output$scorePlot <- renderPlot({
            df<-patient_data
            df$highlight<-"no"
            
            if(isTruthy(input$Patient)){
                df$highlight[df$Patient==input$Patient]<-"yes"
            }
            
            # bargraph 
            x    <- ggplot(data=df, aes(x=Patient, y=Score, fill=df$highlight)) + 
                               geom_bar(stat="identity") +
                               scale_fill_manual(values=c("yes"="tomato", "no"="gray"), guide=FALSE)
            x
            
        })
        
        
        output$scoreTCGA<-renderPlot({
                df<-tcga_data
                df$highlight<-"no"
                
                if(isTruthy(input$TCGA)){
                    df$highlight[df$FileName==input$TCGA]<-"yes"
                }
                
                # bargraph 
                x    <- ggplot(data=df, aes(x=FileName, y=HRDLOH, fill=df$highlight)) + 
                    geom_bar(stat="identity") +
                    scale_fill_manual(values=c("yes"="tomato", "no"="gray"), guide=FALSE)
                x
                
            })
    
    output$mutationCount <- renderPlot({
        # stacked bargraph
        plot    <- ggplot(data=mutation_count, aes(x=mutation_count[,1], y=mutation_count[,2]))+geom_bar(stat="identity")
        plot
    })
    
    output$dt_PatientData<-DT::renderDataTable(server=FALSE, {
        df<-patient_data
        DT::formatRound(DT::datatable(df[, c("Patient","Score")],
                                      extensions = 'Buttons',
                                      
                                      options = list(
                                          paging = TRUE,
                                          searching = TRUE,
                                          fixedColumns = TRUE,
                                          autoWidth = TRUE,
                                          ordering = TRUE,
                                          dom = 'tB',
                                          buttons = c('copy', 'csv', 'excel'),
                                          search = list(
                                              regex = TRUE, 
                                              caseInsensitive = FALSE, 
                                              search = 'M[ae]')
                                      ),
                                      
                                      class = "display", rownames=TRUE),
                       columns=c("Patient", "Score"))
   })
 }

# Run the application 
    shinyApp(ui = ui, server = server)

}

hrd_shiny(patient_data)
