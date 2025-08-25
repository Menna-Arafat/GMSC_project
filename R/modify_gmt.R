
library("BioCor")
library("GSEABase")
library("org.Hs.eg.db")
library(dplyr)

#read gmt
gmt <- geneIds(getGmt("Final Custom copy.gmt" ,
                      geneIdType=SymbolIdentifier()))
head(gmt)

#substitute any space in gene name
gmt_mod= lapply(gmt, function(x) gsub(" +", "-", x))

#export gmt
gmt_table= ldply(gmt_mod, rbind)
gmt_table[is.na(gmt_table)]= ""
write.table(gmt_table, "gmt_mod.gmt", sep = "\t", row.names = FALSE, col.names = FALSE)


#-----------------------------------
#gmt in the long format (a data.frame of 2 column with term and gene)
gmt_long=  gmt_table %>%
  pivot_longer(
    cols= -1,
    names_to = "term",
    values_to = "gene"
  ) %>%
  filter(gene != "" & !is.na(gene)) %>% 
  dplyr::select(-"term")

head(gmt_long$gene, n=100)
write.table(gmt_long, "gmt_long.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

#--------------------------------
#modify structure of enrichment file to get it as # list of pathways,
#names of the list are the pathways names and the elements are the genes
paths= paths %>%  .[which(colnames(.) %in% c("ID", "geneID"))]
paths= paths[!duplicated(paths[,1]), ]
row.names(paths)= NULL
paths= paths %>% column_to_rownames("ID")

#when we have the dataframe with named rows, we can get a named lists by transforming every row to list
paths_lists= apply(paths,1, as.list)
#initialize empty list
paths_lists_splited= vector("list", length= length(paths_lists)) %>% setNames(., names(paths_lists))

for(i in seq_along(paths_lists)){
paths_lists_splited[[i]]  = unlist(lapply(paths_lists[[i]]$geneID, function(s) strsplit(s,"/")))
}


df= stack(paths_lists_splited)
#-----------------------------------