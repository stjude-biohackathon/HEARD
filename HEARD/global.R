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
id_data<- read.table(file.path((paste(directory, "/TCGA_ID_signature_exposures.txt", sep = ''))), sep="\t", header=TRUE)
sbs_data <- read.table(file.path((paste(directory, "/TCGA_SBS_signature_exposure.txt", sep = ''))), sep="\t", header=TRUE)
mut_call_data <- read.table(file.path((paste(directory, "/TCGA_HRD_positive_samples_mut_calls.txt", sep = ''))), sep="\t", header=TRUE)
