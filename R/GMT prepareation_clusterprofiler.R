# ----------------------
# GSEA tutorial
# ----------------------

# Setting up environment ===================================================

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Set seed
set.seed(123456)
install.packages("BiocParallel")

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(fgsea)
library(tidyverse) # includes ggplot2, for data visualization. dplyr, for data manipulation.
library(RColorBrewer) # for a colorful plot
library(pheatmap)
library(clusterProfiler) # for PEA analysis
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot) # for visualizations
library(ggupset) # for visualizations
library(ggtree)
library(ReactomePA)
library(airway)
library(SummarizedExperiment)
library(msigdbr)
library("pathview")
library(qusage)
library(snow)
library(stats)
# Set relevant paths
list.files()
in_path <- "Input/" # input path, where your data is located
out_path <- "Output/" # output path, where you want your results exported to
bg_path <- "GMT/" # folder with your background genes used for PEA


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea")
BiocManager::install('qusage')
BiocManager::install('clusterProfiler')
BiocManager::install('snow')

# Functions ===================================================
## Function: Adjacency matrix to list -------------------------
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

## Function: prepare_gmt --------------------------------------
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  # for debug
  file <- gmt_files[1]
  genes_in_data <- input$ids
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    # Save to GMT format
    output_filename <- paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.gmt')
    write_gmt(final_list, output_filename)
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}
# Analysis ====================================================

## 1. Read in data -----------------------------------------------------------
list.files(in_path)
df <- read.csv(paste0(in_path, 'severevshealthy_degresults.csv'), row.names = 1)

## 2. Prepare background genes -----------------------------------------------

# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
my_genes <- input$ids
list.files(bg_path)
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files
bg_genes <- prepare_gmt(gmt_files[1], my_genes, savefile = TRUE)



# Analysis ====================================================

## 1. Read in data -----------------------------------------------------------
# Get all the genes in your dataset and assign them to my_genes 
input <- read.csv(paste0(in_path, 'Input for testing.csv'))
# Extract the values and names from the data frame
values <- input$gene.list
names <- input$ids
# Create the geneList object with named numeric values
geneList <- as.numeric(values)
names(geneList) <- names
## 2. Prepare background genes -----------------------------------------------
my_genes <- input$ids
list.files(bg_path)

Custom.gmt <- read.gmt(paste0(bg_path, 'custom.gmt'))
str(Custom.gmt)
setwd('C:\\Users\\RSH-DA\\Desktop\\test fgsea\\GMT')
Custom.gmt <- read.csv('term2gene.csv', header=F)
## 4. Run GSEA ---------------------------------------------------------------
# Easy peasy! Run fgsea with the pathways 
RES1 <- GSEA(geneList, TERM2GENE = Custom.gmt, pvalueCutoff = 1 , 
             pAdjustMethod = "fdr", minGSSize = 10 , maxGSSize = 1000,
             exponent = 2)
head(RES1)
# Save as.csv
GSEA.RESULT <- as.data.frame(RES1)
write.csv(GSEA.RESULT, paste0(out_path, 'GSEA-RESULT.csv'))


