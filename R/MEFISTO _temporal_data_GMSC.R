setwd("C:/Users/USER/Documents/Github/GMSC_paper")
dir.create("plots")
dir.create("output")
# BiocManager::install("fields")
library("GSVAdata")
library(MOFA2)   
library("GSVA")
library(GSEABase)
library(GSVAdata)
library(fgsea)
library(tibble)
library("RColorBrewer")
library("circlize")
library(ComplexHeatmap)
library(biomaRt)
library(dplyr)
library(plyr)
library(fields)
library(ggplot2)

# load data
data= read.csv(paste0("input/","Final.DataSet.csv") )
metadata= data.frame(sample= colnames(data)[-1],
                     condition= as.character(data[1,-1]))
data= data[-1,]
head(data)
head(metadata)

# change data type to numeric
data[,-1]= lapply(data[,-1], function(x) as.numeric(as.character(x))) 
sum(is.na(data))

# mapping HMDB id to gene names
mapping= read.delim(paste0("input/","idmapping.tsv" ))
data$id= mapping$To[match(data$Sample, mapping$From)]
length(unique(data$id))

# remove duplicate ids
data= data %>% mutate(mean= rowSums(dplyr::select(.,where(is.numeric)))/ncol(data)-2) %>% 
                group_by(id) %>%
                slice(which.max(mean)) %>% 
                ungroup() %>% select(-c(mean, Sample)) %>% 
                filter(!is.na(id)) %>% as.data.frame()

row.names(data)= data$id
data= data[,!grepl("id", colnames(data))]
head(data)

# set metadata
no_diff= colnames(data)[grepl("7G1", colnames(data))]
early_diff= colnames(data)[grepl("7G3|7G4|7G2", colnames(data))]
late_diff= colnames(data)[grepl("14G", colnames(data))]

metadata$condition %>% unique()
metadata$differentiation.state= ifelse(metadata$sample %in% no_diff, "No differentiation",
                                       ifelse(metadata$sample %in% early_diff, "Early differentiation", "Late differentiation"))

metadata$time.points= ifelse(metadata$sample %in% no_diff, "Timepoint-0 (7G1)",
                             ifelse(metadata$sample %in% early_diff, "Timepoint-1 (7d)", "Timepoint-2 (14d)"))

metadata$time= ifelse(metadata$sample %in% no_diff, 0,
                      ifelse(metadata$sample %in% early_diff, 1, 2))

head(metadata)

# ----------------------------------------------------------------------------------
# MOFA pipeline
data_list= list(as.data.frame(data))
data_list= lapply(data_list, as.matrix)

# create MOFA
MOFAobject <- create_mofa(data_list, groups = NULL, extract_metadata = TRUE)
MOFAobject

# covariate matrix with samples in columns
time = metadata[,c("sample","time"), drop=FALSE] %>% column_to_rownames("sample") %>% t()
rownames(time) <- "time"

# Add sample metadata to the model
samples_metadata(MOFAobject) <- metadata
# set covariates
MOFAobject <- set_covariates(MOFAobject, covariates = time)
MOFAobject

# check data options
data_opts <- get_default_data_options(MOFAobject)
data_opts$center_groups= TRUE # set TRUE if the data is not mean centered
# check model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$ard_factors= TRUE   # ard_factors: use ARD prior in the factors? Default is TRUE if using multiple groups.
# ard_weights: use ARD prior in the weights? Default is TRUE if using multiple views.
model_opts$num_factors= 3
# training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "medium"
train_opts$drop_factor_threshold= -1 #  a value of 0.01 implies that factors explaining less than 1% of variance (in each view) will be dropped
# MEFISTO options
mefisto_opts <- get_default_mefisto_options(MOFAobject)

# prepare object mofa
MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           mefisto_options = mefisto_opts,
                           model_options = model_opts,
                           training_options = train_opts
)

# train mofa
#MOFAobject <- run_mofa(MOFAobject, outfile= "output/MOFA2_GMSC.hdfs", use_basilisk = TRUE)

#saveRDS(MOFAobject,"output/MOFA2_GMSC.rds")
# ------------------------------------------------------------------------------
# read output mofa object
MOFAobject = readRDS("output/MOFA2_GMSC.rds")

# plot varance explained by factors
p= plot_variance_explained(MOFAobject)
print(p)
ggsave("plots/variance_explained_by_factors_heatmap.png", p,  height = 4, width = 10,   dpi = 600)

# scale parameters for each factor, which give us an indication of the smoothness per factor 
#along the covariate (here time) and are between 0 and 1. A scale of 0 means that 
#the factor captures variation independent of time, a value close to 1 tells us that this factor
#varies very smoothly along time.
get_scales(MOFAobject)

plot_factors_vs_cov(MOFAobject, color_by = "time")

df <- plot_factors_vs_cov(MOFAobject, color_by = "time",
                          legend = FALSE, return_data = TRUE)
head(df)
df=df[df$factor== "Factor1",]

# line plot for factor 1
p= ggplot() +
  geom_point(data = df, aes(x = value.covariate, y = value.factor)) +
  geom_smooth(data = df, aes(x = value.covariate, y = value.factor), 
              method = "gam",formula = y ~ s(x, k = 3),  color =  "#FFA500CC" ,fill ="#F3F586FF", linewidth = 1.5)+
  scale_x_continuous(breaks = c(0,1,2),
                     labels = c("Timepoint-0 (7G1)","Timepoint-1 (7d)","Timepoint-2 (14d)"))+
  labs(title = "Covariation of MOFA factor 1 with time points|phases of differentiation", 
       x = "", y = "Latent Factor 1") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA),  
    panel.grid.major = element_line(color = "grey90"),            
    panel.grid.minor = element_line(color = "grey90"),
    plot.title = element_text(size = 11)
  )

print(p)
ggsave("plots/Covariation_MOFA_factor1_with_time_points_gam.png", p,  height = 4, width = 10,   dpi = 600)
# ------------------------------------------------------------------------------
# Top Weight features
p= plot_top_weights(MOFAobject, factors = 1,nfeatures = 30, view = 1)+
  labs(title="Top 30 features for MOFA factor1" ) +
  theme(
   strip.background = element_blank(),
   strip.text = element_blank(),
   plot.title = element_text( size = 13))
print(p)
ggsave("plots/top_30_features_factor1.png", p,  height = 9, width = 8,   dpi = 600) 


# function to extract top features with regards to its weight
get_top_features= function(weight_mtx, direction=NULL, top_n= NULL, factors= NULL){
  
  factor_weights= lapply(factors, function(x) {
    weight = as.data.frame(weight_mtx)
    weight$symbol = rownames(weight)
    
     if(is.null(direction)){
      weight= weight %>% arrange(desc(abs(.[[x]]))) %>% slice(1:top_n) 
    }else if(direction== "pos"){                          
      weight= weight %>% arrange(desc(.[[x]])) %>% slice(1:top_n)
    }else if(direction== "neg"){
      weight= weight %>% arrange(.[[x]]) %>% slice(1:top_n) 
    }
    
    # symbol=  gsub("_Methylation|_mRNA|_Proteome", "",  row.names(weight))
    symbol= weight$symbol
    return(symbol)  })
  
  return(factor_weights)    
}

weight= get_weights(MOFAobject)$view_1
feature.pos.w = get_top_features(weight_mtx= weight,direction="pos", top_n= 30, factors= c(1) ) %>% unlist()
                                                                 

head(feature.pos.w) 
write.csv(feature.pos.w , "output/top_30_features_factor1.csv")
