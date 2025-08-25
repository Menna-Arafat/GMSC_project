#tutoria:https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html#:~:text=Gene%20set%20variation%20analysis%20(GSVA,from%20genes%20to%20gene%20sets.

setwd("C:/Users/USER/Documents/mesenchymal_differentiation/ssgsea")
dir.create("plots")
dir.create("output")

#' ## BiocManager::install("fields")
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
library(limma)
library(fields)
library(tidyr)
#' ## load your gmt as gene set (list of vectors)
list.files("C:/Users/USER/Documents/resources/")
#' gset = gmtPathways("C:/Users/USER/Documents/resources/Human_GO_AllPathways_withPFOCR_no_GO_iea_July_01_2024_symbol.gmt.txt")
#' #' ## filter gene set to have at least 5 terms pr set
#' gset_sub= gset[lapply(gset, length) >= 5]
#' 
#' gset_main= gset_sub[!grepl("PMC",names(gset_sub))]
#' gset_pmc= gset_sub[grepl("PMC",names(gset_sub))]
#' writeGmtPathways(gset_main, "C:/Users/USER/Documents/resources/Human_GO_AllPathways_main.gmt.txt")
#-----------------------------------------------------------------------------------------

gset = gmtPathways("C:/Users/USER/Documents/resources/Human_GO_AllPathways_main.gmt.txt")
#' ## filter gene set to have at least 5 terms pr set
gset_sub= gset[lapply(gset, length) >= 5]
#' ##  convert list to long formats for clusterprofiler
gmt_long= stack(gset_sub) 
gmt_long= gmt_long[,c(2,1)] 


#' ## load data
list.files()
data= read.csv("Final.DataSet.csv") 
metadata= data.frame(sample= colnames(data)[-1],
                     condition= as.character(data[1,-1]))
data= data[-1,]
#' ## change type
data[,-1]= lapply(data[,-1], function(x) as.numeric(as.character(x))) 
sum(is.na(data))
str(data)
#' ## id mapping
list.files()
mapping= read.delim("idmapping.tsv" )
data$id= mapping$To[match(data$Sample, mapping$From)]
length(unique(data$id))
#' ## remove duplicate ids
data= data %>% mutate(mean= rowSums(dplyr::select(.,where(is.numeric)))/ncol(data)-2) %>% 
                         group_by(id) %>%
                         slice(which.max(mean)) %>% 
                         ungroup() %>% dplyr::select(-c(mean, Sample)) %>% 
                         filter(!is.na(id)) %>% as.data.frame()
row.names(data)= data$id
data= data[,!grepl("id", colnames(data))]
head(data)

#' ## determine the distribution of data usually gaussian for log transformed data, and poisson for FPKM count data
general_mean= apply(data, 1, mean)
hist(general_mean, main = "Histogram of Data", xlab = "Values", breaks = 30)
#' ## check normality
shapiro.test(general_mean)
#' ## log transform
data.log= log(data, base=2)
#' ## check normality
general_mean= apply(data.log, 1, mean)
shapiro.test(general_mean)
#' ## check whether it follows gaussian distribution 
qqnorm(general_mean)

#' ## check whether it follows poisson distribution 
general_mean =apply(data.log, 1, mean) #' ## general_mean=lambda
var= apply(data.log, 1, var)
cor(lambda_est, var)
poisson_quantiles = qpois(ppoints(length(general_mean)), lambda = general_mean)  #' ##  Theoretical quantiles, #' ## propability density function of poisson model probability distribution that models the number of times an event happens in a fixed interval of time or space given the average rate is constant, quantile function is essentially the inverse of this CDF, as it finds the k-value corresponding to a specific cumulative probability

#' ## plot the empirical quantiles against the theoretical quantiles.
qqplot(poisson_quantiles, sort(general_mean), 
       main = "Q-Q Plot for Poisson Distribution",
       xlab = "Theoretical Quantiles (Poisson)",
       ylab = "Empirical Quantiles (Data)")


#' ## run ssgsea
gsva_obj = gsvaParam(as.matrix(data),
                      gset_sub, minSize=1, maxSize=500,
                      kcdf="auto")
ES = gsva(gsva_obj)
write.csv(ES, "single_sample_gsea_enrichment_scores.csv", row.names = T)
#---------------------------------------------------------------
#---------------------------------------------------------------
data= read.csv("Final.DataSet.csv") 
metadata= data.frame(sample= colnames(data)[-1],
                     condition= paste0("X", as.character(data[1,-1])))

metadata$condition= factor(metadata$condition, levels = c(paste0("X7G",1:4), paste0("X14G", 1:4) ) )
ES= read.csv("single_sample_gsea_enrichment_scores.csv") %>% column_to_rownames("X")

#' ## Differential expression at pathway level (apply limma on pathway enrichment scores between groups/subgroups)
design= model.matrix(~ 0+metadata$condition) %>% as.data.frame()
row.names(design)= metadata$sample
names(design)= gsub("metadata\\$condition", "", names(design))
colnames(design) <- make.names(colnames(design))
colnames(design) 

library(limma)
fit <- lmFit(ES, design)
fit <- eBayes(fit)
res <- decideTests(fit, p.value=0.01)
summary(res)
#' ## DE pathways
coef= fit$coefficients %>% as.data.frame()
tt <- topTable(fit, coef=4, n=Inf)
DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.05]
DEpwys_es = ES[DEpwys, ]

write.csv(tt, "output/DE_pathways_limma_X7G4.csv", row.names = T)
write.csv(DEpwys_es , "output/DEpwys_ES_X7G4.csv", row.names = T)

#------------------------------
##  Get the results for a specific contrast
coef(fit) %>% head()
contrast <- makeContrasts(X7G1 - X14G4, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

##  Get the top differentially expressed genes
tt <- topTable(fit2, number = Inf)
DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.05]
DEpwys_es = ES[DEpwys, ]

write.csv(tt, "output/DE_pathways_limma_X7G1-X14G4.csv", row.names = T)
write.csv(DEpwys_es , "output/DEpwys_ES_X7G1-X14G4.csv", row.names = T)


#--------------------------------------------------
#' ## GSVA scores have higher precision for larger gene sets
gssizes <- geneSetSizes(ES)
plot(sqrt(gssizes), sqrt(fit$sigma), xlab="Sqrt(gene sets sizes)",
     ylab="Sqrt(standard deviation)", las=1, pch=".", cex=3)
#' ## In such a setting, we can improve the analysis of differentially expressed pathways by using the limma-trend approach (Phipson et al. 2016)
#' ## setting the trend parameter in the call to the eBayes() function to the vector of gene set sizes. 
#' ##  Key Idea Behind Limma-Trend
#' ##  In RNA-seq and pathway analysis, there can be a mean-variance trend—i.e.,
#' ##  genes or pathways with higher expression tend to have larger variances.
#' ##  The limma-trend approach takes this into account by allowing the prior variance estimate
#' ##  to depend on the mean expression of the gene or pathway, rather than assuming a constant variance for all genes.
fit <- eBayes(fit, trend=gssizes)
res <- decideTests(fit, p.value=0.01)
summary(res)


#' ## -------------------------------------------------------------------------------
#' ## heatmap DE pathways
DEpwys_es=read.csv( "output/DEpwys_ES_X7G4_selected.csv" ) %>% column_to_rownames("X")
DEpwys_es= DEpwys_es[!grepl("LEISHMANIA|PARASITIC|GOCC|GOMF", row.names(DEpwys_es)), ] 


#' ## for hierarchical clustering,to determine whether to use pearson (assume norrmally distributed data and linear relation) or spearman (assume non normally distributed data)
shapiro.test(as.numeric(ES[1, ]))  #' ## p-value <= 0.05 indicative of non normal distribution 
shapiro.test(ES[, 1]) 

colorLegend <- c("#FFA500","#8C510A" , "darkslateblue","darkolivegreen" ,"#F6E8C3" ,"#DFC27D" , "steelblue","#ADFF2F" )
names(colorLegend) <- unique(metadata$condition)
sample.color.map = colorLegend[metadata$condition]
names(sample.color.map) <- metadata$sample

sampleClustering <- hclust(as.dist(1-cor(as.matrix(DEpwys_es), method="spearman")), #' ## pearson #' ## spearman
                           method="complete")
geneSetClustering <- hclust(as.dist(1-cor(t(DEpwys_es), method="spearman")),
                            method="complete")
palette <- colorRampPalette(c( "lightyellow2" ,"#DFC27D" , "#FCAA0FFF", "darkred"))(256)



png("plots/heatmap_DE_pathsways_xxxX7G4.png", height=4100, width = 2650, res= 600)
heatmap(as.matrix(DEpwys_es), ColSideColors=sample.color.map, xlab="samples",
        ylab="", margins=c(2, 20),
        col = palette ,
        labRow=substr(gsub("_", " ", gsub("^KEGG_|^REACTOME_|^BIOCARTA_|\\%GO:.*", "", rownames(DEpwys_es))), 1, 35),
        labCol="", 
        scale="none", 
        Colv= as.dendrogram(sampleClustering),
        Rowv=as.dendrogram(geneSetClustering)
         )
legend("topright", inset = 0.001, names(colorLegend), fill=colorLegend, bg="white",
       cex = .6, title= "Conditions") #' ##  #' ## x = .01, y =9
image.plot(zlim = range(DEpwys_es, na.rm = TRUE),
           legend.only = TRUE, 
           horizontal = TRUE, 
           legend.shrink = 0.3, 
           legend.width = 0.9,
           legend.mar= 4.5,
           col = palette , 
           legend.position = c(0.5, 0.5), 
           inset = c(-1.1, 0),
           legend.args = list(text = "Enrichment Score (ES)", side = 3, line = .5, cex = 0.6))
dev.off()
#' ## -------------------------------------------------------------------------------
#' ## volcano plot for DE pathways
tt <- topTable(fit, coef=2, n=Inf)
DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.01]


plot(tt$logFC, -log10(tt$P.Value), pch=16, cex=1, col=grey(0.75),
     main="", xlab="GSVA enrichment score difference", las=1,
     ylab=expression(-log[10]~~Raw~P-value))
abline(h=-log10(max(tt$P.Value[tt$adj.P.Val <= 0.01])),
       col=grey(0.5), lwd=1, lty=2)
points(tt$logFC[match(DEpwys, rownames(tt))],
       -log10(tt$P.Value[match(DEpwys, rownames(tt))]),
       pch=16, cex=1, col="darkred")
text(max(tt$logFC)*0.85, -log10(max(tt$P.Value[tt$adj.P.Val <= 0.01])),
     "1% FDR", pos=3)
#' ## -------------------------------------------------------------------------------
#' ## heatmap enrichment scores 
heatmap_data <- scale(as.matrix(ES))

Heatmap(
  matrix = as.matrix(heatmap_data ),
  name = "Enrichment Scores",
  col = colorRamp2(c(-3, 0, 3), c("#9370DB", "#F7F5F4","orange" )),#' ## matlab::jet.colors(200),
  show_row_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = T
  #' ## top_annotation  = ta
  #' ## column_title = ""
) 




