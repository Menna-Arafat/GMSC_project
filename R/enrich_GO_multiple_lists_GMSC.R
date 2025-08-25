setwd(r'{C:\Users\USER\Documents\GMSC\go_ontology}')


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
#get and prepare gmt 
# gmt= list()
# for (file in list.files(pattern= ".gmt")){
#   x = geneIds(getGmt(file,  geneIdType=SymbolIdentifier()))
#   gmt[[file]]= x
# }
 
# #concatnate lists
# combined_list <- do.call(c, gmt)
#filter sets with less than 5 terms per set 
# gmt_sub <- combined_list[lengths(combined_list) > 4]
# gmt_sub %>% head()

gmt= gmtPathways("c:/Users/USER/Documents/resources/Human_GOBP_AllPathways_noPFOCR_with_GO_iea_July_01_2024_symbol.gmt.txt")
gmt_go = geneIds(getGmt(r'{C:\Users\USER\Documents\resources\c5.go.v2024.1.Hs.symbols.gmt}' , geneIdType=SymbolIdentifier()))

gmt= c(gmt, gmt_go)
gmt_long= stack(gmt) 
gmt_long= gmt_long[,c(2,1)]  

#convert gmt to long format so that every term has its corresponding gene
# gmt_long=  gmt %>%
#             pivot_longer(
#               cols= -1,
#               names_to = "term",
#               values_to = "gene"
#             ) %>%
#             filter(gene != "" & !is.na(gene)) %>% 
#             dplyr::select(-"term")
#-------------------------------------------------------------------------------
#load data
id_list= list()
for (file in list.files(pattern= ".csv")){
  x= read.csv(file)
  name= gsub("\\s+\\(1\\)\\.csv", "", file)
  proteins= x$X[x$p.value <= .05] %>% unlist() %>% unname()
  id_list[[name]]= proteins 
}

#convert id
mapping_file=  read.delim("idmapping.tsv")


#mapping function
map_id= function ( vector, mapping_file){
  mapping= as.data.frame(mapping_file)
  names(mapping)= c("query", "name")
  mapped_vec= mapping$name[match(vector , mapping$query)] %>% .[!is.na(.)] %>% unique()
  return(mapped_vec)
}

mapped_list <- lapply(id_list, map_id, mapping_file = mapping_file)
#-------------------------------------------------------------------------------
enrichment= "ORA"  # "GSEA" , "GO"

if(enrichment== "ORA"){
#enrichment multiple lists (ORA)
enrich= llply(mapped_list, function(i) {
          x= enricher(i, TERM2GENE =gmt_long, pvalueCutoff = 1, minGSSize = 5, qvalueCutoff = 1)
          return(as.data.frame(x))})

}else if(enrichment== "GSEA" ){
#enrichment multiple lists (GSEA)
enrich= llply(mapped_list, function(i) {
  x= GSEA(FC_vector, TERM2GENE = gmt_long, pvalueCutoff = 0.05,qvalueCutoff = 1, minGSSize = 5)
  return(as.data.frame(x))})

}else if(enrichment== "GO" ){
#enrichment multiple lists(GO)
enrich= llply(mapped_list, function(i) {
  x= enrichGO(gene =i,
              universe = NULL, #background genes
              OrgDb ='org.Hs.eg.db', 
              keyType = "SYMBOL", #SYMBOL, #ENTREZID 'UNIPROT'
              readable = T,
              ont = "ALL",
              minGSSize = 5,
              pvalueCutoff =.5, 
              qvalueCutoff = 1)
  return(as.data.frame(x))
  })
}

enrich.df = ldply(enrich, rbind)
enrich.df= enrich.df[enrich.df$pvalue <= .05, ]
enrich.df$ID= gsub("%.*", "", enrich.df$ID)
enrich.df$Description= gsub("%.*", "", enrich.df$Description)
enrich.df <- plyr::rename(enrich.df, c(.id = "Cluster"))
enrich.df$Cluster = as.factor(enrich.df$Cluster)
write.csv(enrich.df, "output/enrichment_all_conditions_clusterprofiler.csv", row.names = F)
#-------------------------------------------------------------------------------
# enrichplot
enrich.df= read.csv("output/enrichment_all_conditions_clusterprofiler.csv")
enrich.df$Cluster = as.factor(enrich.df$Cluster)
#let the dataframe output inherit the class of compareclusterResult/ clusterprofiler
res <- new("compareClusterResult", compareClusterResult = enrich.df, 
           geneClusters = mapped_list , .call = match.call(expand.dots = TRUE))

# Then create the color vector
clusters <- unique(res@compareClusterResult$Cluster)
col <- palt(length(clusters))
names(col) <- clusters

# Check consistency
print(names(col))
print(clusters)

#visualize enrichmap
palt = colorRampPalette(c( "#BA6756", "#B0C4DE", "#E3B31C",
                          "tan", "purple",  "maroon3" , "#21908CFF" , "grey66"))
col= palt(length(unique(enrich.df$Cluster)))
names(col)= unique(enrich.df$Cluster)

p= cnetplot(res, max.overlaps= 50) +
  theme_minimal()+
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),            
    axis.line = element_blank(),            
    axis.text = element_blank(),            
    axis.ticks = element_blank(),             
    axis.title = element_blank()
  )+
  labs(fill= "Groups") +
  scale_color_manual(values = col)

print(p)
ggsave("plots/enrich_map_wgcna_modules.png", p, width= 28, height=19)

#---------------------------------------