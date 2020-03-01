#step 1: read data####
path_machine = "/Users/dcongtinh/gene-abundance/experiment/results/fc_model_with_feature_selected/qtf_pc576_10_fillseqf_nb10_auy_gray/"
path_r =path_machine
setwd(path_r)
dataset_name <- "T2dgene"
file_name <- paste0('feature_selected_', tolower(dataset_name), '.csv')
table = read.csv(file_name, header = TRUE)
dim(table)
#step 2: visualize feature selected####
library(ggplot2)
title_size <- 30
font_size <- 20
ggplot(data=table, aes(x=freq)) + ylab("Count")+ xlab("Frequency of selected times") + 
   geom_bar() + stat_count(aes(y=..count..,label=..count..),geom="text",vjust=-0.8,size=5.5) +
      theme(axis.text=element_text(size=font_size),
          axis.title=element_text(size=font_size,face="bold"),
          plot.title = element_text(hjust = 0.5, size=title_size, face="bold", margin=margin(8, 0, 16, 0))) +
              ggtitle(paste0("Histogram for selected ", dataset_name, "'s features"))
