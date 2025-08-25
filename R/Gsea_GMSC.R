setwd(r'{C:\Users\USER\Documents\mesenchymal_differentiation\GSEA}')

dir.create("input")
dir.create("output")
dir.create("plots")
#load libraries
library(plyr)
library(dplyr)
library(igraph)
library(tidyr)
library(readxl)
library(openxlsx)
library(org.Hs.eg.db)
library(scales)
library(ggplot2)
library(tibble)
library("GSEABase")
library("EnrichmentBrowser")
library(clusterProfiler)
library(enrichplot)
library(biomaRt)
library(fgsea)

gmt= gmtPathways("c:/Users/USER/Documents/resources/Human_GOBP_AllPathways_noPFOCR_with_GO_iea_July_01_2024_symbol.gmt.txt")

gmt_long= stack(gmt) 
gmt_long= gmt_long[,c(2,1)]  

#' ## load data
list.files()
data= read.csv("input/Final.DataSet.csv") 
metadata= data.frame(sample= colnames(data)[-1],
                     condition= as.character(data[1,-1]))
data= data[-1,]
#' ## change type
data[,-1]= lapply(data[,-1], function(x) as.numeric(as.character(x))) 
sum(is.na(data))
str(data)
#-------------------------------------------------------------------------------
#' ## id mapping
list.files()
mapping= read.delim("input/idmapping.tsv" )
data$id= mapping$To[match(data$Sample, mapping$From)]
length(unique(data$id))
#' ## remove duplicate ids
data= data %>% mutate(mean= rowSums(dplyr::select(.,where(is.numeric)))/ncol(data)-2) %>% 
                group_by(id) %>%
                slice(which.max(mean)) %>% 
                ungroup() %>% dplyr::select(-c(mean, Sample)) %>% 
                filter(!is.na(id)) %>% as.data.frame()

row.names(data)= NULL
data= data %>% column_to_rownames("id")
head(data)
#-------------------------------------------------------------------------------
#prepare FC
groups = list(c("7G1", "7G2"),
             c("7G1" , "7G3"),
             c("7G1", "7G4" ),
             c("14G1", "14G2"),
             c("14G1" , "14G3"),
             c("14G1", "14G4" ),
             c("7G2", "14G2"),
             c("7G3" , "14G3"),
             c("7G4", "14G4"))

group_names <- sapply(groups, function(x) paste(x, collapse = "_"))

get_FC= function(data, group, direction= NULL, FC= NULL){

  data= as.data.frame(data)
  FC= as.numeric(as.character(FC))
  G1= group[1]
  G2= group[2]
  group1_means = rowMeans(data[,grepl(G1, colnames(data))], na.rm = T)
  group2_means = rowMeans(data[,grepl(G2, colnames(data))], na.rm = T)

  #let the denominator 1 to avoid division by zero
  fold_change = ifelse(group1_means == 0, group2_means / (group1_means + 1), group2_means / group1_means)
  # 0/0 returns NA, so Set fold_change to 0 where both group means are zero
  fold_change[is.na(fold_change)] = 0
  names(fold_change)=  gsub("_Methylation|_mRNA|_Proteome|_lnRNA|_miRNA", "",  names(fold_change))
  if(is.null(direction)){
    fold_change= fold_change
  }else if(direction== "UP"){
    fold_change= fold_change[fold_change > FC]
  }else if( direction== "DOWN"){
    fold_change= fold_change[fold_change < FC]
  }

  return(fold_change)
}


fc <- lapply(groups, function(group) {
  get_FC(data = data, group = group, direction = NULL, FC = NULL)
})
names(fc) <- group_names

#-------------------------------------------------------------------------------
#enrich clusterprofiler
enrich=  lapply(fc, function(i){
  #GSEA
  i=  i[order(i, decreasing = TRUE)]
  i = i[!duplicated(names(i))]  
  e = GSEA(i, TERM2GENE = gmt_long, pvalueCutoff = 1, minGSSize = 1)
  #ORA
  #e= enricher(i, TERM2GENE =TF, minGSSize = 5, pvalueCutoff = 1,  qvalueCutoff = 1)
  return(as.data.frame(e))
})

sig= lapply(enrich, function(i){
  i= i[i$pvalue <= .1, ]
  return(as.data.frame(i))
})

sig_df=ldply(sig, rbind)
sig_df= sig_df[!grepl("%GO|LEISHMANIA|PARASIT*", sig_df$Description),]

#write.csv(sig_df, "output/gsea_all_groups.csv", row.names = F)

sig= read.csv( "output/gsea_all_groups.csv")
sig= sig[sig$pvalue <= .1, ]
write.csv(sig, "output/gsea_all_groups_.1pvalue.csv", row.names = F)

