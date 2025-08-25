#' ---
#' title: "Pathway Crosstalk Analysis"
#' author: "Menna Arafat"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---

#' ## Set Working Directory and Create Directories
setwd("C:/Users/USER/Documents/mesenchymal_differentiation")
dir.create("plots")


#' ## Load Required Libraries
library(plyr)
library(dplyr)
library(igraph)
library(tidyr)
library(readxl)
library(openxlsx)
library(scales)
library(ggplot2)
library(tibble)
library(future)
library(future.apply)
library(biomaRt)
library(enrichplot)
library(ggrepel)
library(ggraph)
library(tidygraph)
library(RColorBrewer)

#' ## Load Data
#' Load the crosstalk data and terms from Excel files.
data = read.xlsx("data/crosstalk_results_permutation.xlsx", sheet = "crosstalk")
terms = read.xlsx("data/crosstalk_results_permutation.xlsx", sheet = "terms") %>%
  column_to_rownames(".id") %>%
  t() %>%
  as.data.frame()

#' ## Filter Significant Pathways
#' Apply a p-value cutoff of 0.1 to select significant pathways.
data = data %>% filter(p.value <= 0.1 & pathway2 != "Disease")
terms = terms[, colnames(terms) %in% c(data$pathway1, data$pathway2)]

#' ## Add Transcription Factors
TF = read_excel("output/enrichment_top_30_features_factor1_gprofiler.xlsx", sheet= "TF")[c(2,11)]
TF$term_name = TF$term_name %>% lapply(., function(i) strsplit(i, ":|;")[[1]][2])

TF_list = TF %>% 
  column_to_rownames("term_name") %>%
  t() %>%
  as.data.frame() %>%
  as.list() %>%
  lapply(., function(x) strsplit(x, ",")[[1]])

#' ## Build Lists of Significant Pathways
paths_lists = as.list(terms)
paths_lists = lapply(paths_lists, function(x) x[x != "" & !is.na(x)])
all_lists = c(paths_lists, TF_list)
all_genes = unlist(paths_lists) %>% unname()

#' ## Convert List to Long Format
list2df <- function(inputList) {
  ldf <- lapply(seq_len(length(inputList)), function(i) {
    data.frame(categoryID = rep(names(inputList[i]), length(inputList[[i]])),
               Gene = inputList[[i]])
  })
  ldply(ldf, rbind)
}

#' ## Convert List to Graph Object
list2graph <- function(inputList, directed = FALSE) {
  x <- list2df(inputList)
  g <- graph_from_data_frame(x, directed = directed)
  return(g)
}

#' ## Build Igraph Object
graph = list2graph(all_lists)

#' ## Set Colors for Igraph
palette = colorRampPalette(c("#F6E8C3", "#B0C4DE", "#E3B31C", "tan", "purple", "maroon3", "#21908CFF", "darkgrey"))
palette(40)
colors = palette(length(all_lists))

all_col = rep("lightgrey", length(V(graph)))

for(i in seq_along(paths_lists)) {
  all_col[V(graph)$name %in% names(paths_lists)[i]] = colors[i]
}
all_col[V(graph)$name %in% names(TF_list)] = "#96D73F"

#' ## Create Tidygraph Object
graph_tbl <- as_tbl_graph(graph) %>%
  activate(nodes) %>%
  mutate(size = ifelse(name %in% names(paths_lists), 6,
                       ifelse(name %in% names(TF_list), 4, 2)),
         label_size = ifelse(name %in% names(paths_lists), 2.5, 2),
         shape = ifelse(name %in% names(TF_list), "2", "1"),
         col = all_col)

#' ## Calculate Centroids for Each Community
layout <- layout_with_kk(graph_tbl)

centroids <- sapply(seq_along(paths_lists), function(i) {
  nodes = V(graph)$name %in% names(paths_lists[[i]])
  colMeans(layout[nodes, ])
}) %>% t() %>% as.data.frame()

row.names(centroids) = names(paths_lists)

#' ## Plot Using Ggraph
p = ggraph(graph_tbl, layout = 'circle') +  
  geom_edge_fan(aes(alpha = 0.5), color = "lightgrey", show.legend = FALSE) + 
  geom_node_point(aes(size = size, color = col, shape = shape), show.legend = FALSE) + 
  scale_color_identity() + 
  scale_size_identity() +
  geom_node_text(aes(label = name, size = label_size), show.legend = FALSE, repel = TRUE) +
  theme_graph(fg_text_colour = 'black') +
  theme(element_text(size = 13))

#' ## Save the Plot
ggsave("plots/crossTalk_circular.png", p, width = 11, height = 7, dpi = 600)


