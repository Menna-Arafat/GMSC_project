setwd("C:/Users/USER/Documents/Github/GMSC_paper")
dir.create("plots")
dir.create("output")

#load libraries
#devtools::install_github("Moonerss/CIBERSORT")
library("RColorBrewer")
library("circlize")
library(dplyr)
library(plyr)
library(CIBERSORT)
library(readxl)
library(tibble)
library(tidyverse)
library(ggplot2)

#load your bulk cell data expression
data= read.csv("input/Final.DataSet.csv") 
metadata= data.frame(sample= colnames(data)[-1],
                     condition= as.character(data[1,-1]))
data= data[-1,]
#change type
data[,-1]= lapply(data[,-1], function(x) as.numeric(as.character(x))) 
sum(is.na(data))

#id mapping
mapping= read.delim("input/idmapping.tsv" )
data$id= mapping$To[match(data$Sample, mapping$From)]
length(unique(data$id))


#remove duplicate ids
data= data %>% mutate(mean= rowSums(dplyr::select(.,where(is.numeric)))/ncol(data)-2) %>% 
                group_by(id) %>%
                slice(which.max(mean)) %>% 
                ungroup() %>% 
                select(-c(mean, Sample)) %>% 
                filter(!is.na(id)) %>% as.data.frame()

head(data)

# Prepare data input for Cibersort
matrix_file= cbind(data$id, data[,-ncol(data)]) %>% as.data.frame()
names(matrix_file)[1] = "Gene.symbol"
str(matrix_file)
write.table(matrix_file, "output/matrix_file.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#load signature matrix (average gene expression per cluster_ single cell reference)
sig_matrix = read_excel("input/sig_matrix_mouse incisor atlas.xls") %>% as.data.frame()
names(sig_matrix)[1]= "Gene.symbol"
str(sig_matrix)
write.table(sig_matrix, "output/sig_matrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#------------------------------------------------------------------------------
#RUN Cibersort algorithm input provided as path to file :))
results <- cibersort(sig_matrix= "output/sig_matrix.txt"   ,
                     mixture_file= "output/matrix_file.txt" )


head(results)
write.csv(results, "output/cibersort_output.csv")

# Visualize
# Convert data into long format for ggplot
# results= read.csv( "output/cibersort_output.csv") %>% column_to_rownames("X")
results= results[, -c(15,16,17)]
results$samples= row.names(results)
res= rbind(results[13:24,], results[1:12,])
data_long= res %>% pivot_longer(names_to = "cell_type",  values_to = "proportion", -samples)
#to get samples ordered on x axis
data_long$samples = factor(data_long$samples, levels = unique(data_long$samples))
head(data_long)  


colors= c("#440154FF" ,"#3B528BFF","maroon2", "#21908CFF", "#DAA520" ,
          "#5DC863FF", "#FDE725FF", "purple", "#DFC27D" ,"#F6E8C3", "#C7EAE5", 
          "#F4A460", "grey", "orange","darkgrey")  
names(colors)= colnames(res)

p= ggplot(data_long, aes(x = samples, y = proportion, fill = cell_type)) +
  scale_fill_manual( values= colors) +
  geom_bar(stat = "identity") +
  labs(x = "Samples", y = "Proportion", fill = "Cell Type") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12), 
    axis.text.y = element_text( size = 15), 
    panel.background = element_rect(fill = "white"),  
    plot.background = element_rect(fill = "white"),   
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank()
    
  ) 
print(p)
ggsave("plots/Cibersort_single_cell_deconvolution.png",p, dpi = 600, width = 11, height = 13)
