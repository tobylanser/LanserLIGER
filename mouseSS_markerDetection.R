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

Sys.setenv('R_MAX_VSIZE'=memory.size(max=TRUE))

seurat3x <- readRDS("seurat3x.rds")

tissue.avgs <- AverageExpression(seurat3x, features = "Foxd1", group.by ='Tissue')
write.csv(tissue.avgs, "seurat3x_Foxd1_tissueExpression.csv")


# Marker detection
seurat3x_cluster_markers <- FindAllMarkers(seurat3x)
seurat3x_cluster_markers_top10 <- seurat3x_cluster_markers %>%
  group_by(cluster) %>% 
  top_n(10, avg_log2FC)
write.csv(seurat3x_cluster_markers_top10, "seurat3x_cluster_markers_top10_dims20_res.4.csv")



# Set the identity of your cells to the desired column

Idents(object = seurat3x) <- seurat3x@meta.data$"orig.ident"

seurat3x_orig.ident_markers <- FindAllMarkers(seurat3x)
seurat3x_orig.ident_markers_top10 <- seurat3x_orig.ident_markers %>%
  group_by(cluster) %>% 
  top_n(50, avg_log2FC)
write.csv(seurat3x_orig.ident_markers_top10, "seurat3x_orig.ident_markers_top50_dims20_res.4.csv")

# Set the identity of your cells to the desired column
Idents(object = seurat3x) <- seurat3x@meta.data$"Tissue"

seurat3x_tissue_markers <- FindAllMarkers(seurat3x)
seurat3x_tissue_markers_top10 <- seurat3x_tissue_markers %>%
  group_by(cluster) %>% 
  top_n(50, avg_log2FC)
write.csv(seurat3x_tissue_markers_top10, "seurat3x_tissue_markers_top50_dims20_res.4.csv")

seurat3xMarkers <- seurat3x
# set cutoff to zero expression
foxd1.pos <- subset(seurat3xMarkers, subset = Foxd1 > 0.32)
foxd1.neg <- subset(seurat3xMarkers, subset = Foxd1 <= 0.32)
WhichCells(seurat3xMarkers@assays[["RNA"]]@data["Foxd1",])>1
subset(seurat3xMarkers, subset = Foxd1 > 0.32)
#create metadata foxd1 column
fodx1.pos.cells <- "Foxd1+"
seurat3x@meta.data <- cbind(seurat3x@meta.data, fodx1.pos.cells)
colnames(seurat3x@meta.data)[which(names(seurat3x@meta.data) == "fodx1.pos.cells")] <- "foxd1"

table(seurat3x@meta.data$foxd1)
seurat3xMarkers <- merge(foxd1.pos, y = foxd1.neg, project = "Foxd1Markers")
seurat3x@meta.data$foxd1 <- seurat3xMarkers@meta.data$foxd1

fodx1.neg.cells <- "Foxd1-"
foxd1.neg@meta.data <- cbind(foxd1.neg@meta.data, fodx1.neg.cells)
colnames(foxd1.neg@meta.data)[which(names(foxd1.neg@meta.data) == "fodx1.neg.cells")] <- "foxd1"

seurat3x <- readRDS("seurat3x.rds")

# Set the identity of your cells to the desired column
Idents(object = seurat3x) <- seurat3x@meta.data$"foxd1"

seurat3x_tissue_markers <- FindAllMarkers(seurat3x)

seurat3x_foxd1_markers <- seurat3x_tissue_markers
seurat3x_foxd1_markers_top10 <- seurat3x_foxd1_markers %>%
  group_by(cluster) %>% 
  top_n(10, avg_log2FC)
write.csv(seurat3x_foxd1_markers_top10, "seurat3x_foxd1_markers_top10_dims20_res.4.csv")

# Pathway analysis
head(AverageExpression(seurat3x, features = "foxd1"))

saveRDS(seurat3x, "seurat3x.rds")

