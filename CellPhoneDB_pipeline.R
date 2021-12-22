##Load CRAN Packages
library(patchwork)
library(cluster)
library(cowplot)
library(plot3D)
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
library(sigora)
library(colorpatch)
library(Seurat)
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
library(viridis)
library(qs)
library(scCustomize)
library(Cairo)
library(iTALK)
library(nichenetr)
library(garnett)
library(Hmisc)
library(survival)
library(lattice)

setwd("/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/rds_files")

adult_ILC_Tcells_fibros <- qread("Adult_ILC_Tcells_fibros_timecourse.rds")


DimPlot_scCustom(adult_ILC_Tcells_fibros)

ct <- as.matrix(adult_ILC_Tcells_fibros$cell_type)

colnames(ct) <- c("cell_type")

write.table(ct, '/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/adult_ILC_Tcells_fibros_celltype.txt', 
            sep = '\t', row.names = T, col.names = T, quote = F)

df <- as.matrix(GetAssayData(object = adult_ILC_Tcells_fibros, slot = "counts"))

rownames(df) <- toupper(rownames(df))

write.table(df, '/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/GeneXCell_adult_ILC_Tcells_fibros.txt', 
            sep = '\t', row.names = T, col.names = T, quote = F)



###Subsetting
b1_b2_all_nonimmune <- qread("timecourse_all_nonimmune.rds", label = T)

Idents(b1_b2_all_nonimmune) <- "seurat_clusters"
FeaturePlot_scCustom(b1_b2_all_nonimmune, features = c("Loxl1", "Col1a1"), label = T, split.by = "Age")
aged_b1_b2_all_nonimmune <- subset(b1_b2_all_nonimmune, idents = "Aged")

FeaturePlot_scCustom(b1_b2_all_nonimmune, features = "Cxcl12", label = T)


DimPlot_scCustom(b1_b2_all_nonimmune) + DarkTheme()

ggsave("title_slide.jpeg", device = "jpeg",
       width = 10, height = 10, units = "in", dpi = 300)

aged_b1_b2_all_nonimmune <- NormalizeData(object = aged_b1_b2_all_nonimmune)
aged_b1_b2_all_nonimmune <- FindVariableFeatures(object = aged_b1_b2_all_nonimmune)
aged_b1_b2_all_nonimmune <- ScaleData(object = aged_b1_b2_all_nonimmune)

###PCA
aged_b1_b2_all_nonimmune <- RunPCA(object = aged_b1_b2_all_nonimmune)
print(x = aged_b1_b2_all_nonimmune[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = aged_b1_b2_all_nonimmune, dims = 1:6)

DimPlot(object = aged_b1_b2_all_nonimmune)
aged_b1_b2_all_nonimmune <- ProjectDim(object = aged_b1_b2_all_nonimmune)

##JackStraw
aged_b1_b2_all_nonimmune <- JackStraw(object = aged_b1_b2_all_nonimmune, num.replicate = 100)
aged_b1_b2_all_nonimmune <- ScoreJackStraw(object = aged_b1_b2_all_nonimmune, dims = 1:20)
JackStrawPlot(object = aged_b1_b2_all_nonimmune, dims = 1:20)

ElbowPlot(object = aged_b1_b2_all_nonimmune, ndims = 20)

###
aged_b1_b2_all_nonimmune <- FindNeighbors(object = aged_b1_b2_all_nonimmune, dims = 1:8)
aged_b1_b2_all_nonimmune <- FindClusters(object = aged_b1_b2_all_nonimmune, resolution = 0.7)
##
aged_b1_b2_all_nonimmune <- RunUMAP(aged_b1_b2_all_nonimmune, dims = 1:8)

FeaturePlot(aged_b1_b2_all_nonimmune, features = "Loxl1", label = T, order = T)

Idents(aged_b1_b2_all_nonimmune) <- "seurat_clusters"
aged_fibros <- subset(aged_b1_b2_all_nonimmune, idents = c("1", "2", "3", "11", "14"))

aged_fibros <- NormalizeData(object = aged_fibros)
aged_fibros <- FindVariableFeatures(object = aged_fibros)
aged_fibros <- ScaleData(object = aged_fibros)

###PCA
aged_fibros <- RunPCA(object = aged_fibros)
print(x = aged_fibros[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = aged_fibros, dims = 1:6)

DimPlot(object = aged_fibros)
aged_fibros <- ProjectDim(object = aged_fibros)

##JackStraw
aged_fibros <- JackStraw(object = aged_fibros, num.replicate = 100)
aged_fibros <- ScoreJackStraw(object = aged_fibros, dims = 1:20)
JackStrawPlot(object = aged_fibros, dims = 1:20)

ElbowPlot(object = aged_fibros, ndims = 20)

###
aged_fibros <- FindNeighbors(object = aged_fibros, dims = 1:5)
aged_fibros <- FindClusters(object = aged_fibros, resolution = 0.7)
##
aged_fibros <- RunUMAP(aged_fibros, dims = 1:5)

FeaturePlot(aged_fibros, features = "Col1a1", label = T, order = T)

aged_fibros$cell_type <- "Fibros"
qsave(aged_fibros, "/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/rds_files/Aged_fibros_timecourse.rds")

