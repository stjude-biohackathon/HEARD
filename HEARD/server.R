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
  
  
      output$LOH_Info<-renderPlot({
        # Pre-processing dataset
        # create a helper variable for sorting
        test <- gene_loh %>% group_by(Gene) %>% count()
        test$sum <- test$n
        test$n <- NULL
        # group by gene and LOH event
        gene_loh_dat <- gene_loh %>% group_by(Gene,LOH_Event) %>% count()
        # convert to numeric for sorting
        gene_loh_dat$Count <- as.numeric(gene_loh_dat$n)
        gene_loh_dat$value <- as.numeric(gene_loh_dat$n)
        # sort by decreasing value
        order_data <- gene_loh_dat[order(-gene_loh_dat$value),]
        # set Gene as a factor for sorting
        order_data$Gene <- as.factor(order_data$Gene)
        dat <- order_data %>% select(Gene, LOH_Event, value)
        dat <- as.data.frame(dat)
        dat <- dat[order(dat$value,decreasing=TRUE),]
        dat_sum <- merge(dat,test,by="Gene")
        
        # Code for the plot without TP53
        #Not in
        `%!in%` = Negate(`%in%`)
        p <- ggplot(dat_sum[dat_sum$Gene %!in% "TP53",],aes(x=reorder(Gene,-sum),y=value,fill=LOH_Event)) +
          geom_bar(stat="identity") +
          scale_fill_brewer(palette = "Paired") +
          theme(
            axis.text.x = element_text(angle = 90,vjust=1, face = "bold"),
            #axis.text.x = element_text(face = "bold"),
            panel.background = element_rect(fill = "transparent"), # bg of the panel
            plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
            panel.grid.major = element_blank(), # get rid of major grid
            panel.grid.minor = element_blank(), # get rid of minor grid
            legend.background = element_rect(fill = "transparent"), # get rid of legend bg
            legend.box.background = element_rect(fill = "transparent"))#, # get rid of legend panel bg
        p
      })
  
      output$id_htmap<-renderPlot({
        mat<-id_score[,2:19]
        rownames(mat)<-id_score[,1]
        chosen <-input$TCGA
        row.colors <- rep.int(c("black"), nrow(mat))
        row.colors[rownames(mat) == chosen] <- "red"
        col.colors <- rep.int(c("black"), ncol(mat))
        col.colors[colnames(mat) == "ID6"] <- "red"
        rownames(mat)<-substr(rownames(mat),1,12)
        col<-colorRampPalette(brewer.pal(8, name="YlGnBu"))(25)
        
        Heatmap(mat, name="ID Exposure", 
                width = ncol(mat)*unit(3, "mm"), 
                height = nrow(mat)*unit(6, "mm"),
                cluster_rows=FALSE, 
                cluster_columns=FALSE,
                column_names_gp = grid::gpar(fontsize = 8, col=col.colors),
                row_names_gp = grid::gpar(fontsize = 8, col=row.colors), col=col)
        
      })
      
      output$sbs_htmap<-renderPlot({
        mat<-sbs_score[,2:79]
        rownames(mat)<-sbs_score[,1]
        chosen <-input$TCGA
        row.colors <- rep.int(c("black"), nrow(mat))
        row.colors[rownames(mat) == chosen] <- "red"
        col.colors <- rep.int(c("black"), ncol(mat))
        col.colors[colnames(mat) %in% c("SBS3","SBS5")] <- "red"
        rownames(mat)<-substr(rownames(mat),1,12)
        col<-colorRampPalette(brewer.pal(8, name="YlGnBu"))(25)
        
        #Heatmap
        Heatmap(mat, name="SBS Exposure", 
                width = ncol(mat)*unit(2, "mm"), 
                height = nrow(mat)*unit(5, "mm"),
                cluster_rows=FALSE, 
                cluster_columns=FALSE,
                column_names_gp = grid::gpar(fontsize = 6, col=col.colors),
                row_names_gp = grid::gpar(fontsize = 8, col=row.colors), col=col)
        
      })
      
      
      output$scoreTCGA<-renderPlot({
        df<-tcga_data
        df$highlight<-"no"
        
        if(isTruthy(input$TCGA)){
          df$highlight[df$FileName==input$TCGA]<-"yes"
        }
        df$FileName<-substr(df$FileName, 1, 12)
        
        x    <- ggplot(data=df, aes_string(x=paste0("reorder(FileName, -", input$Column,")"), 
                                           y=input$Column, fill="highlight")) + 
          geom_bar(stat="identity") +
          theme(axis.text.x = element_text(angle = 90))+
          scale_fill_manual(values=c("yes"="tomato", "no"="gray"), guide=FALSE)
        x
        
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
