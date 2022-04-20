library("shiny")
library("tidyverse")
library("readxl")
library("readr")
library("dplyr")
library("tools")
library("DT")
#library("bslib")
library("shinycustomloader")
library("plotly")
library("shinyWidgets")
library("ggplot2")
library("RColorBrewer")
suppressPackageStartupMessages(library("ComplexHeatmap"))
library("R.utils")
library("ncdf4")
options(shiny.maxRequestSize=300*1024^2)
options(repos = BiocManager::repositories())

#load example dataset first
directory <- (file.path("example_dataset"))
#SMD read-in
tcga_data<-read.table(file.path((paste(directory, "/TCGA_HRD_positive_samples.txt", sep = ''))), sep="\t", header=TRUE)
sbs_score<-read.table(file.path((paste(directory, "/TCGA_SBS_signature_exposure.txt", sep = ''))), sep="\t", header=TRUE)
id_score<-read.table(file.path((paste(directory, "/TCGA_ID_signature_exposures.txt", sep = ''))), sep="\t", header=TRUE)
gene_loh<-read.table(file.path((paste(directory, "/gene_LOH_events.txt", sep = ''))), sep="\t", header=TRUE)

#EVAN read-in
clinical_data <- read.table(file.path((paste(directory, "/TCGA_HRD_positive_samples_mut_calls.txt", sep = ''))), sep="\t", header=TRUE)
hrd_data <- read.table(file.path((paste(directory, "/TCGA_HRD_positive_samples.txt", sep = ''))), sep="\t", header=TRUE)
id_data <- read.table(file.path((paste(directory, "/TCGA_ID_signature_exposures.txt", sep = ''))), sep="\t", header=TRUE)
sbs_data <- read.table(file.path((paste(directory, "/TCGA_SBS_signature_exposure.txt", sep = ''))), sep="\t", header=TRUE)
mut_call_data <- read.table(file.path((paste(directory, "/TCGA_HRD_positive_samples_mut_calls.txt", sep = ''))), sep="\t", header=TRUE)


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



function(input, output, session) {
  
##### Dan's code below##################################
  #data table upload functions moved to top for clarity (replace example dataset with uploaded tar.gz)
  uploaded_HRD_data <- eventReactive(c(input$uploaded_HRD_data),
                                     { if (file_ext(input$uploaded_HRD_data$datapath) == "gz") {
                                       isGzipped(input$uploaded_HRD_data$datapath)
                                       tarFile <- gunzip(input$uploaded_HRD_data$datapath)
                                       unlink("temp", recursive = T, force = T)
                                       untar(tarFile, exdir = "temp/")
                                       directory <- list.dirs("temp", recursive = FALSE)
                                       
                                       #SMD data read-in
                                       tcga_data<<-read.table(file.path((paste(directory, "/TCGA_HRD_positive_samples.txt", sep = ''))), sep="\t", header=TRUE)
                                       updatePickerInput(session = session, inputId = "TCGA",
                                                         choices = (tcga_data[,1]))
                                       updatePickerInput(session = session, inputId = "Column", selected = "HRDsum",
                                                         choices=colnames(tcga_data)[!colnames(tcga_data) %in% c("FileName")])
                                       sbs_score<<-read.table(file.path((paste(directory, "/TCGA_SBS_signature_exposure.txt", sep = ''))), sep="\t", header=TRUE)
                                       id_score<<-read.table(file.path((paste(directory, "/TCGA_ID_signature_exposures.txt", sep = ''))), sep="\t", header=TRUE)
                                       gene_loh<<-read.table(file.path((paste(directory, "/gene_LOH_events.txt", sep = ''))), sep="\t", header=TRUE)

                                       #EVAN read-in
                                       clinical_data <<- read.table(file.path((paste(directory, "/TCGA_HRD_positive_samples_mut_calls.txt", sep = ''))), sep="\t", header=TRUE)
                                       hrd_data <<- read.table(file.path((paste(directory, "/TCGA_HRD_positive_samples.txt", sep = ''))), sep="\t", header=TRUE)
                                       updatePickerInput(session = session, inputId = "hrd_metrics", selected = "HRDsum",
                                                         choices=colnames(hrd_data[,-1]))
                                       id_data <<- read.table(file.path((paste(directory, "/TCGA_ID_signature_exposures.txt", sep = ''))), sep="\t", header=TRUE)
                                       sbs_data <<- read.table(file.path((paste(directory, "/TCGA_SBS_signature_exposure.txt", sep = ''))), sep="\t", header=TRUE)
                                       updatePickerInput(session = session, inputId = "pat_id",
                                                         choices=(sbs_data$FileName))
                                       mut_call_data <<- read.table(file.path((paste(directory, "/TCGA_HRD_positive_samples_mut_calls.txt", sep = ''))), sep="\t", header=TRUE)
                                       
                                       
                                       list_of_df <- list(hrd_data, id_data, sbs_data)
                                       
                                       data_viz <- Reduce(function(x, y) merge(x, y, all=TRUE), list_of_df, accumulate=FALSE)
                                       
                                       data_viz
                                     }
                                     })
##########################################################################
  # Evan function to render input dynamically
#  output$data <- renderTable({
#    file <- input$file1
#    ext <- tools::file_ext(file$datapath)
    #
#    req(file)
#    validate(need(ext=="txt","Please upload a txt file"))
#    tcga_data <- read.delim(file$datapath, sep="\t", header=TRUE)
#  })
  #observe({print(input$pat_id)})
  # patient data text renderer
  output$clinical_info <- renderUI({
    #clinical_data$Sample_type
    df <- clinical_data[clinical_data$Patient == substr(input$pat_id, 1, 12),]
    if (length(df$Sample_type) > 1) {
      # check the sample type
      print(df$Sample_type)
      # split rows into two data frames
      df1 <- df[1,]
      df2 <- df[2,]
      #check if first row contains info for germline
      # we need to turn this into a function later
      if(df1$Sample_type == "Normal") {
        print(paste("first sample type is normal:",df1$Sample_type))
        text0 <- paste("",sep="<br/>")
        text1 <- paste("Patient Name:",df1$Patient)
        text2 <- paste("Patient Diagnosis:",df1$Dx)
        text20 <- paste("",sep="<br/>")
        text3 <- paste("Patient Germline Mutation Information:")
        text4 <- paste("Gene: ",df1$annovar_region_gene)
        text5 <- paste("Chromosome: ",df1$Chr)
        text6 <- paste("Position: ",df1$Pos)
        text7 <- paste("Type: ",df1$Type)
        text8 <- paste("dbSNP:",df1$dbSNP)
        text21 <- paste("% Alt. Allele:", df1$Percent_alternative_allele)
        text9 <- paste("Class:",df1$Class)
        text10 <- paste("Function:",df1$annovar_exonic_function)
        #HTML(paste(text0,text1,text2,text3,text4,sep="<br/>"))
      }
      #check if first row contains info for tumor
      else if (df1$Sample_type == "Tumor") {
        print(paste("first sample type is Tumor:",df1$Sample_type))
        text0 <- paste("",sep="<br/>")
        text1 <- paste("Patient Name:",df1$Patient)
        text2 <- paste("Patient Diagnosis:",df1$Dx)
        text20 <- paste("",sep="<br/>")
        text3 <- paste("Patient Somatic Mutation Information:")
        text4 <- paste("Gene: ",df1$annovar_region_gene)
        text5 <- paste("Chromosome: ",df1$Chr)
        text6 <- paste("Position: ",df1$Pos)
        text7 <- paste("Type: ",df1$Type)
        text8 <- paste("dbSNP:",df1$dbSNP)
        text21 <- paste("% Alt. Allele:", df1$Percent_alternative_allele)
        text9 <- paste("Class:",df1$Class)
        text10 <- paste("Function:",df1$annovar_exonic_function)
        #HTML(paste(text0,text1,text2,text3,text4,sep="<br/>"))
      }
      #check if second row contains info for germline
      if(df2$Sample_type == "Normal") {
        text11 <- paste("",sep="<br/>")
        text12 <- paste("Patient Germline Mutation Information:")
        text13 <- paste("Gene: ",df2$annovar_region_gene)
        text14 <- paste("Chromosome: ",df2$Chr)
        text15 <- paste("Position: ",df2$Pos)
        text16 <- paste("Type: ",df2$Type)
        text17 <- paste("dbSNP:",df2$dbSNP)
        text22 <- paste("% Alt. Allele:", df2$Percent_alternative_allele)
        text18 <- paste("Class:",df2$Class)
        text19 <- paste("Function:",df2$annovar_exonic_function)
        #HTML(paste(text0,text1,text2,text3,text4,sep="<br/>"))
      } 
      #check if second row contains info for tumor
      else if (df2$Sample_type == "Tumor") {
        text11 <- paste("",sep="<br/>")
        text12 <- paste("Patient Somatic Mutation Information:")
        text13 <- paste("Gene: ",df2$annovar_region_gene)
        text14 <- paste("Chromosome: ",df2$Chr)
        text15 <- paste("Position: ",df2$Pos)
        text16 <- paste("Type: ",df2$Type)
        text17 <- paste("dbSNP:",df2$dbSNP)
        text22 <- paste("% Alt. Allele:", df2$Percent_alternative_allele)
        text18 <- paste("Class:",df2$Class)
        text19 <- paste("Function:",df2$annovar_exonic_function)
        #HTML(paste(text0,text1,text2,text3,text4,sep="<br/>"))
      }
      HTML(paste(text0,text1,text2,text20,text3,text4,text5,text6,text7,text8,text21,text9,text10,
                 text11,text12,text13,text14,text15,text16,text17,text22,text18,text19,sep="<br/>"))
    }
    else if (df$Sample_type == "Normal") {
      df1 <- df
      print(paste("sample type is normal:",df1$Sample_type))
      text0 <- paste("",sep="<br/>")
      text1 <- paste("Patient Name:",df1$Patient)
      text2 <- paste("Patient Diagnosis:",df1$Dx)
      text20 <- paste("",sep="<br/>")
      text3 <- paste("Patient Germline Mutation Information:")
      text4 <- paste("Gene: ",df1$annovar_region_gene)
      text5 <- paste("Chromosome: ",df1$Chr)
      text6 <- paste("Position: ",df1$Pos)
      text7 <- paste("Type: ",df1$Type)
      text8 <- paste("dbSNP:",df1$dbSNP)
      text21 <- paste("% Alt. Allele:", df1$Percent_alternative_allele)
      text9 <- paste("Class:",df1$Class)
      text10 <- paste("Function:",df1$annovar_exonic_function)
      HTML(paste(text0,text1,text2,text20,text3,text4,text5,text6,text7,text8,text21,text9,text10,sep="<br/>"))
    }
    else if (df$Sample_type == "Tumor") {
      df1 <- df
      print(paste("first sample type is Tumor:",df1$Sample_type))
      text0 <- paste("",sep="<br/>")
      text1 <- paste("Patient Name:",df1$Patient)
      text2 <- paste("Patient Diagnosis:",df1$Dx)
      text20 <- paste("",sep="<br/>")
      text3 <- paste("Patient Somatic Mutation Information:")
      text4 <- paste("Gene: ",df1$annovar_region_gene)
      text5 <- paste("Chromosome: ",df1$Chr)
      text6 <- paste("Position: ",df1$Pos)
      text7 <- paste("Type: ",df1$Type)
      text8 <- paste("dbSNP:",df1$dbSNP)
      text21 <- paste("% Alt. Allele:", df1$Percent_alternative_allele)
      text9 <- paste("Class:",df1$Class)
      text10 <- paste("Function:",df1$annovar_exonic_function)
      HTML(paste(text0,text1,text2,text20,text3,text4,text5,text6,text7,text8,text21,text9,text10,sep="<br/>"))
    }
  })
  # output chomosome images
  output$chromImage <- renderImage(deleteFile = FALSE,{
    filename <- normalizePath(file.path(paste('example_dataset/TCGA_HRD_positive_samples_CNA_figs/',
                                              input$pat_id,'omosome_view/',
                                              input$pat_id,'omosome_view',
                                              chr_translation[input$Chromosome],
                                              '.png', sep='')))
    #print(filename)
    list(src = filename,height=800,width=1000)
  })
  # Function for rendering hrdscores
  output$hrdScores<-renderPlot({
    # set df to our hrd metrics
    df <- hrd_data
    # get averages
    avg <- dplyr::summarize_all(df[,-1],mean)
    avg$FileName <- "Average"
    df_plot <- rbind(df[df$FileName == input$pat_id,],
                     avg)
    # set up highlight
    df_plot$Sample <- "no"
    df_plot$Sample <- ifelse(df_plot$FileName == input$pat_id,"yes","no")
    df_plot$FileName <- substr(df_plot$FileName, 1, 12)
    metric <- input$hrd_metrics
    print(df_plot$Sample)
    #ggplot(hrd_data, aes(x=FileName,y=input$hrd_metrics)) +
    ggplot(df_plot, aes(x=FileName,fill=Sample)) +
      #y=input$hrd_metrics, fill="highlight")) +
      geom_bar(stat="identity",aes_string(y=input$hrd_metrics))+
      #coord_flip() +
      xlab("Patient") +
      ylab("Value") +
      #theme_bw(base_size = 16) +
      scale_fill_manual(values=c("yes"="tomato", "no"="gray"))+
      theme(
        axis.text.x = element_text(angle = 90,vjust=1, face = "bold",size=12),
        axis.text.y = element_text(face = "bold",size=12),
        axis.title = element_text(face = "bold",size=14),
        legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(face="bold",size=12),
        #axis.text.x = element_text(face = "bold"),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))
  })
  # function for mutsigs plot
  output$mutsigs <- renderPlot({
    mat <- cbind(sbs_data[,2:79],
                 id_data[,2:19])
    row.names(mat) <- sbs_data[,1]
    chosen <- input$pat_id
    mat <- mat %>%
      dplyr::select(SBS3,SBS5,ID6)
    #row.names(mat) <- mat$FileName
    row.colors <- rep.int(c("black"), nrow(mat))
    row.colors[rownames(mat) == chosen] <- "red"
    col.colors <- rep.int(c("black"), ncol(mat))
    col.colors[colnames(mat) %in% c("SBS3","SBS5","ID6")] <- "red"
    rownames(mat)<-substr(rownames(mat),1,12)
    col<-colorRampPalette(brewer.pal(8, name="YlGnBu"))(25)
    mat <- mat[row.names(mat) == getShortName(input$pat_id),]
    #Heatmap
    Heatmap(as.matrix(mat), name="Exposure", 
            width = ncol(mat)*unit(3, "cm"), 
            height = nrow(mat)*unit(4, "cm"),
            cluster_rows=FALSE, 
            cluster_columns=FALSE,
            column_names_gp = grid::gpar(fontsize = 14, col=col.colors),
            column_title_gp = grid::gpar(fontsize = 16, col=col.colors),
            row_names_gp = grid::gpar(fontsize = 14, col=row.colors), col=col,
            row_title_gp = grid::gpar(fontsize = 16, col=col.colors),
            show_heatmap_legend = FALSE,
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.1f", mat[i, j]), x, y, gp = gpar(fontsize = 16,col="tomato"))})
    
  })
  #### Evan Code Here####################################################################################
  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      row.names(hrd_data) <- hrd_data$FileName
      row.names(sbs_data) <- sbs_data$FileName
      row.names(id_data) <- id_data$FileName
      # remove duplicate columns
      sbs_data$FileName <- NULL
      id_data$FileName <- NULL
      # create binded table
      params <- cbind(hrd_data[row.names(hrd_data) == as.character(input$pat_id),],
                      sbs_data[row.names(sbs_data) == as.character(input$pat_id),],
                      id_data[row.names(id_data) == as.character(input$pat_id),])
      #
      # reduce params to those of interest
      params <- params %>%
        dplyr::select(FileName,SBS3,SBS5,ID6,HRDsum,cellularity)
      # create params list variable
      params_list <- list(
        FileName=params$FileName,
        HRDsum=params$HRDsum,
        cellularity=params$cellularity,
        SBS3=params$SBS3,
        SBS5=params$SBS5,
        ID6=params$ID6,
        HRDStatus=ifelse(input$hrd_stat == "Yes","HRD Positive","HRD Negative"))
      #
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params_list,
                        envir = new.env(parent = globalenv())
      )
    }
  )

  
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
          scale_fill_manual(values=c("yes"="tomato", "no"="gray"), guide="none")
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
