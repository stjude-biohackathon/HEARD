# Gene vs LOH barchart

# Dataset
gene_loh <- read.table(file = "gene_LOH_events.txt",header = TRUE)
#
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
# Code fo the plot with all genes:
p <- ggplot(dat_sum,aes(x=reorder(Gene,-sum),y=value,fill=LOH_Event)) +
  geom_bar(stat="identity") +
  scale_fill_brewer(palette = "Paired") +
  theme(
    axis.text.x = element_text(angle = 90, face = "bold"),
    #axis.text.x = element_text(face = "bold"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"))#, # get rid of legend panel bg
#
# Code for the plot without TP53
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
#
ggsave(p, filename = "gene_vs_loh_barplot_sorted_noTP53.pdf", width=15,bg = "transparent",device="pdf")
#