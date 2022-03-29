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
library(cli)
library(sctransform)
library(scCustomize)
library(CellChat)

setwd("C:/Users/Stevens Lab/Desktop/LanserLIGER/data")

Sys.setenv('R_MAX_VSIZE'=memory.size(max=TRUE))

all_ages.cells <- readRDS("timecourse_allAges_allCells.rds")

Idents(all_ages.cells) <- "Age"

adult.p5 <- subset(all_ages.cells, idents = "Aged", invert = T)

adult.p5 <- NormalizeData(adult.p5)
adult.p5 <- FindVariableFeatures(adult.p5)
adult.p5 <- ScaleData(adult.p5)
###PCA
adult.p5 <- RunPCA(object = adult.p5, verbose = F)
ElbowPlot(object = adult.p5, ndims = 50)
###
adult.p5 <- FindNeighbors(object = adult.p5, dims = 1:40)
adult.p5 <- FindClusters(object = adult.p5, resolution = 0.8)
adult.p5 <- RunUMAP(adult.p5, dims = 1:40)

DimPlot_scCustom(adult.p5, label = T, raster = F, group.by = "cell_type")


saveRDS(adult.p5, "adult_p5_all_cells.rds")

Idents(adult.p5) <- "cell_type"
adult.p5_macs.fibros <- subset(adult.p5, idents = c("Macrophages", "Fibroblasts"))
adult.p5_macs.fibros <- NormalizeData(adult.p5_macs.fibros)
adult.p5_macs.fibros <- FindVariableFeatures(adult.p5_macs.fibros)
adult.p5_macs.fibros <- ScaleData(adult.p5_macs.fibros)
###PCA
adult.p5_macs.fibros <- RunPCA(object = adult.p5_macs.fibros, verbose = F)
ElbowPlot(object = adult.p5_macs.fibros, ndims = 50)
###
adult.p5_macs.fibros <- FindNeighbors(object = adult.p5_macs.fibros, dims = 1:40)
adult.p5_macs.fibros <- FindClusters(object = adult.p5_macs.fibros, resolution = 0.8)
adult.p5_macs.fibros <- RunUMAP(adult.p5_macs.fibros, dims = 1:40)

DimPlot_scCustom(adult.p5_macs.fibros, label = T, raster = F, group.by = "Age")

saveRDS(adult.p5_macs.fibros, "adult.p5_macs.fibros.rds")

adult.p5_macs <- subset(adult.p5, idents = c("Macrophages"))
adult.p5_macs <- NormalizeData(adult.p5_macs)
adult.p5_macs <- FindVariableFeatures(adult.p5_macs)
adult.p5_macs <- ScaleData(adult.p5_macs)
###PCA
adult.p5_macs <- RunPCA(object = adult.p5_macs, verbose = F)
ElbowPlot(object = adult.p5_macs, ndims = 50)
###
adult.p5_macs <- FindNeighbors(object = adult.p5_macs, dims = 1:40)
adult.p5_macs <- FindClusters(object = adult.p5_macs, resolution = 0.8)
adult.p5_macs <- RunUMAP(adult.p5_macs, dims = 1:40)

DimPlot_scCustom(all_cell, label = T, raster = F, group.by = "cell_type", split.by = "Age")

Idents(all_cell) <- "cell_type"
FeaturePlot(adult.p5_fibros, features = c("Ramp2", "Ackr3"), split.by = "Age")

saveRDS(adult.p5_macs, "adult.p5_macs.rds")

adult.p5_fibros <- subset(adult.p5, idents = c("Fibroblasts"))
adult.p5_fibros <- NormalizeData(adult.p5_fibros)
adult.p5_fibros <- FindVariableFeatures(adult.p5_fibros)
adult.p5_fibros <- ScaleData(adult.p5_fibros)
###PCA
adult.p5_fibros <- RunPCA(object = adult.p5_fibros, verbose = F)
ElbowPlot(object = adult.p5_fibros, ndims = 50)
###
adult.p5_fibros <- FindNeighbors(object = adult.p5_fibros, dims = 1:40)
adult.p5_fibros <- FindClusters(object = adult.p5_fibros, resolution = 0.8)
adult.p5_fibros <- RunUMAP(adult.p5_fibros, dims = 1:40)

Idents(adult.p5) <- "seurat_clusters"
DimPlot_scCustom(adult.p5_fibros, label = T, raster = F)

saveRDS(adult.p5_fibros, "adult.p5_fibros.rds")
