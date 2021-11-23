#writing in Robj files from tabula muris
#un-comment LIGER package line if correct

library(Seurat)
library(tidyseurat)
library(rliger)
library(tidyverse)
library(reticulate)
library(ggbeeswarm)
library(ggrepel)
library(ggridges)
library(RColorBrewer)
library(devtools)
library(ggplot2)
library(dplyr)
library(Matrix)
library(umap)
library(cluster)
library(cowplot)
library(Rtsne)
library(scatterplot3d)
library(plotly)
library(ggsci)
library(ggrepel)
library(ggiraph)
library(gplots)
library(ggridges)
library(ggstance)
library(autoplotly)
library(RColorBrewer)
library(pheatmap)
library(treemap)
library(FateID)
library(colorpatch)
library(patchwork)


setwd("C:/Users/Stevens Lab/Desktop/LanserLIGER")

Sys.setenv('R_MAX_VSIZE'=8e+11)
Sys.getenv('R_MAX_VSIZE')

DuraCalcStrom <- readRDS("Dura_calv_stromal_2.rds")
StromFibMyo <- readRDS("StromFibroMyo_clustered_TissueRenamed.rds")
mouseSSfibro <- readRDS("Mouse_SS_Fibro.RDS")

mouseSSfibro$orig.ident <- "Turley"
StromFibMyo$orig.ident <- "Tabula Muris"
DuraCalcStrom$orig.ident <- "Meninges"

fibro.big <- readRDS("fibroMerged.rds")

#rename project.name
DuraCalcStrom@project.name <- "DuraCalcStrom"
StromFibMyo@project.name <- "TabulaMuris"
mouseSSfibro@project.name <- "FibroblastSet"

#Rename tissue slot to Tissue
#StromFibMyo@meta.data$Tissue <- StromFibMyo@meta.data$tissue
#rm(StromFibMyo@meta.data$tissue)


DuraCalcStrom <- UpdateSeuratObject(object = DuraCalcStrom)
mouseSSfibro <- UpdateSeuratObject(object = mouseSSfibro)
StromFibMyo <- UpdateSeuratObject(object = StromFibMyo)

#saveRDS(StromFibMyo, "StromFibroMyo_clustered_TissueRenamed.rds")


rm(DuraCalcStrom)
rm(mouseSSfibro)
rm(StromFibMyo)
rm(Liger3X)
rm(seurat3x)
rm(seurat3xScaled)


fibro.big <- merge(DuraCalcStrom, y = c(mouseSSfibro, StromFibMyo), add.cell.ids = c("DuraCalcStrom", "Shannon Turley", "TabulaMuris"), project = "MergedFibro")


fibro.big <- ScaleData(fibro.big)
fibro.big <- FindVariableFeatures(fibro.big)
fibro.big <-RunPCA(fibro.big, dims = 1:5)

fibro.big <-RunUMAP(fibro.big, dims = 1:5)
fibro.big <-RunTSNE(fibro.big, check_duplicates = FALSE)

#saveRDS(fibro.big, "fibroMerged.rds") 

fibro.big <- readRDS("fibroMerged.rds")

fibro.big@raw.data
DuraCalcStrom@raw.data
Liger3xtake2 <- seuratToLiger(fibro.big, combined.seurat = TRUE, assays.use = "RNA", meta.var = "orig.ident")



Liger3xtake2 <- seuratToLiger(list(DuraCalcStrom, StromFibMyo, mouseSSfibro), names =c("DuraCalcStrom","TabulaMuris", "Shannon Turley"))



head(Liger3xtake2@cell.data)
head(DuraCalcStrom@meta.data)

#mouseSSFibro annotations
SSfibro_df=as.data.frame(mouseSSfibro@meta.data)
row.names(SSfibro_df)=gsub("1_1","",row.names(SSfibro_df))


SSfibro_df=as.data.frame(mouseSSfibro@meta.data)
row.names(SSfibro_df)=gsub("1_1$","",row.names(SSfibro_df))
row.names(SSfibro_df)=gsub("_1_1$","",row.names(SSfibro_df))
head(SSfibro_df)

#DuraCalcStrom annotations
DuraCalcStrom_df=as.data.frame(DuraCalcStrom@meta.data)
row.names(DuraCalcStrom_df)=gsub("1_1","",row.names(DuraCalcStrom_df))


DuraCalcStrom_df=as.data.frame(DuraCalcStrom_df@meta.data)
row.names(DuraCalcStrom_df)=gsub("1_1$","",row.names(DuraCalcStrom_df))
row.names(DuraCalcStrom_df)=gsub("_1_1$","",row.names(DuraCalcStrom_df))
head(DuraCalcStrom_df)

#Tabula Muris annotations
StromFibMyo_df=as.data.frame(StromFibMyo@meta.data)
row.names(StromFibMyo_df)=gsub("1_1","",row.names(StromFibMyo_df))


StromFibMyo_df=as.data.frame(StromFibMyo_df@meta.data)
row.names(StromFibMyo_df)=gsub("1_1$","",row.names(StromFibMyo_df))
row.names(StromFibMyo_df)=gsub("_1_1$","",row.names(StromFibMyo_df))
head(StromFibMyo_df)


Liger3X <-  readRDS("3xLIGER_preNormalizeRenamed.rds")
Liger3X <- normalize(Liger3X)
Liger3X <- selectGenes(Liger3X)
Liger3X <- scaleNotCenter(Liger3X)


#Merged seurat dataset annotations
mergedFibro_df=as.data.frame(fibro.big@meta.data)
row.names(mergedFibro_df)=gsub("1_1","",row.names(mergedFibro_df))


mergedFibro_df=as.data.frame(fibro.big@meta.data)
row.names(mergedFibro_df)=gsub("1_1$","",row.names(mergedFibro_df))
row.names(mergedFibro_df)=gsub("_1_1$","",row.names(mergedFibro_df))
head(mergedFibro_df)

attributes(Liger3xtake2)@cell.data$meta.data <- mergedFibro_df

# Liger3X <-  readRDS("3xLIGER_preNormalizeRenamed.rds")
Liger3xtake2 <- normalize(Liger3xtake2)
Liger3xtake2 <- selectGenes(Liger3xtake2)
Liger3xtake2 <- scaleNotCenter(Liger3xtake2)

### Stage II: Joint Matrix Factorization
#any k value between 20-40 is appropriate. Higher k should be used w/ datasets w/ more substructure
# k.suggest <- suggestK(strom_liger, num.cores = 1, gen.new = F, plot.log2 = T,
# nrep = 2)

#k.suggest <- suggestK(Liger3X, num.cores = 5, nrep = 5)


Liger3xtake2 <- optimizeALS(Liger3xtake2, k = 20)

# saveRDS(Liger3xtake2, "Liger3xtake2_postOptimizeALS.rds")
# 
# Liger3X <- readRDS("3xLIGER_postOptimizeALS.rds")

### Stage III: Quantile Normalization and Joint Clustering
#The `quantile_norm` procedure produces joint clustering assignments and a low-dimensional representation
#that integrates the datasets together. These joint clusters directly from inmf can be used for downstream
#analyses (see below). Alternatively, you can also run Louvain community detection, an algorithm commonly
#used for single-cell data, on the normalized cell factors. The Louvain algorithm excels at merging small clusters
#into broad cell classes and thus may be more desirable in some cases than the maximum factor assignments produced directly by inmf.

### Stage III: Quantile Normalization and Joint Clustering

Liger3xtake2 <- quantile_norm(Liger3xtake2)

Liger3xtake2 <- louvainCluster(Liger3xtake2, resolution = 1)

### Stage IV: Visualization and Downstream Analysis

Liger3xtake2 <- runUMAP(Liger3X)

Liger3X <- readRDS("3xLIGER_postUMAP.rds")

seurat3x <- ligerToSeurat(Liger3X)

seurat3x <- ScaleData(seurat3x)


seurat3x <- FindNeighbors(seurat3x, reduction = "inmf", dims = 1:ncol(Embeddings(seurat3x, "inmf"))) %>%
  FindClusters(resolution = 3)
seurat3x <- RunUMAP(seurat3x, reduction = "inmf", dims = 1:ncol(seurat3x[["inmf"]]))

saveRDS(seurat3x, "seurat3x.rds")








saveRDS(Liger3xtake2, "3xLIGER_postUMAP.rds")

Liger3X <- readRDS("3xLIGER_postUMAP.rds")



saveRDS(seurat3x, "seurat3x.rds")

seurat3x <-readRDS("seurat3xScaled.rds")

Liger3X$umap <- seurat3x@reductions$umap
seurat3x@meta.data
seurat3x$Tissue <- fibro.big$Tissue

seurat3x$Tissue <-md.tissue
rm(md.tissue)

Liger3X <-readRDS("3xLIGER_postUMAP.rds")



seurat3x@meta.data <- fibro.big@meta.data
seurat3x$orig.ident <- fibro.big$orig.ident



seurat3x_df=as.data.frame(seurat3x_df@meta.data)
row.names(seurat3x_df)=gsub("1_1$","",row.names(seurat3x_df))
row.names(seurat3x_df)=gsub("_1_1$","",row.names(seurat3x_df))
head(mergedFibro_df)

write.table(seurat3x$orig.ident, "seurat3x2_orig.ident.txt")

Liger3xtake2$Tissue <- fibro.big$Tissue





#`plotByDatasetAndCluster` returns two graphs, generated by t-SNE or UMAP in the previous step.
#The first colors cells by dataset of origin, and the second by cluster as determined by Liger.
#The plots provide visual confirmation that the datasets are well aligned and the clusters are consistent
#with the shape of the data as revealed by UMAP.

all.plots <- plotByDatasetAndCluster(Liger3X, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T)
all.plots[[1]] + all.plots[[2]]

p_a <- plotByDatasetAndCluster(Liger3X, return.plots = T) 
p_a[[1]] <- p_a[[1]] + theme_classic() + theme(legend.position = c(0.85, 0.15)) + 
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
print(p_a[[1]])

#To directly study the impact of factors on the clustering and determine what genes load most
#highly on each factor, we use the plotGeneLoadings function, which returns plots of factor loading
#on the dimensionally reduced graphs and highly loaded genes by dataset for each factor.
gene_loadings <- plotGeneLoadings(Liger3X, do.spec.plot = FALSE, return.plots = TRUE)
gene_loadings[[4]]


#Using the `runWilcoxon` function, we can next identify gene markers for all clusters.
#We can also compare expression within each cluster across datasets.
#The function returns a table of data that allows us to determine the significance of each gene's differential expression,
#including log fold change, area under the curve and p-value.
cluster.results <- runWilcoxon(Liger3X, compare.method = "clusters")
head(cluster.results)

datasets.results <- runWilcoxon(Liger3X, compare.method = "datasets")
head(datasets.results)

write.csv(datasets.results, "DatasetsResults.csv")

cluster.results <- cluster.results[cluster.results$padj < 0.05,]
cluster.results <- cluster.results[cluster.results$logFC > 3,]

#You can then re-sort the markers by its padj value in ascending order and choose the top 100 for each cell type.
#For example, we can subset and re-sort the output for Cluster 3 and take the top 20 markers by typing these commands:
wilcoxon.cluster_3 <- cluster.results[cluster.results$group == 3, ]
wilcoxon.cluster_3 <- wilcoxon.cluster_3[order(wilcoxon.cluster_3$padj), ]
markers <- wilcoxon.cluster_3[1:20, ]
head(markers)

#We can then visualize the expression profiles of individual genes, such as the differentially expressed genes that we just identified.
#This allows us to visually confirm the cluster- or dataset-specific expression patterns of marker genes.
#`plotGene` returns graphs of gene loading on the dimensionally reduced graph for each dataset.
PRF1 <- plotGene(Liger3X, "Col1a1", axis.labels = c('tSNE 1', 'tSNE 2'), return.plots = T)
PRF1[[1]]
PRF1[[2]]


#We can also use `plotGene` to compare the loading of cluster markers within and between datasets.
IFIT3 <- plotGene(Liger3X, "Gsta4", axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE)
IFITM3 <- plotGene(Liger3X, "Ly6c1", axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE)
plot_grid(IFIT3[[1]],IFIT3[[2]],IFITM3[[1]],IFITM3[[2]], ncol=2)

#We can also identify some shared and dataset-specific markers for each factor and plot them to help in cluster annotation.
#This information can aid in finding cell-type specific dataset differences. We can return these in table format with `getFactorMarkers`,
#though an easy way to visualize these markers (and distribution of factor loadings) is with the package's `plotWordClouds` function.
### CHANGE - Here we can notice that in the plot of factor 9, which seems to correspond to the NK cell cluster, one of the most highly-loading dataset-specific markers is CD16 (FCGR3A). 

# This function uses the V loadings to identify dataset-specific genes and a combination of the W 
# and V loadings to identify shared markers
#CHANGE DATASET
markers <- getFactorMarkers(Liger3X, dataset1='seurat3xProject', dataset2='Steady_State', num.genes = 10)
# The first three elements in this list are the dataset1-specific, shared, and dataset2-specific 
# markers respectively
# Let's take a look at factor 9
head(markers$Steady_State[markers$Steady_State$factor_num == 9, ])

# Show top 10 highly loading dataset specific and shared genes
# Prevent text progress bar from printing
word_clouds <- plotWordClouds(Liger3X, num.genes = 10, do.spec.plot = F, return.plots = T)
print(word_clouds[[9]])

# Show top 10 highly loading dataset specific and shared genes
# with plots of gene loading values
gene_loadings <- plotGeneLoadings(Liger3X, num.genes = 10, do.spec.plot = F, return.plots = T)
print(gene_loadings[[9]])

#We can now plot this gene to better observe finer-grained dataset specific differences. 
p_g2 <- plotGene(Liger3X, 'Gsta4', return.plots = T)
plot_grid(plotlist = p_g2)


