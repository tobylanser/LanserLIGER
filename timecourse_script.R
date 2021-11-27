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


setwd("/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data")

b1_b2_all <- readRDS("Meninges_B1_B2_filtered_updated.rds")

b1_b2_all <- UpdateSeuratObject(object = b1_b2_all)

saveRDS(b1_b2_all, "/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/Meninges_B1_B2_filtered_updated.rds")

DimPlot(b1_b2_all, reduction = "umap", pt.size = 0.75, label = TRUE, label.size = 5)
FeaturePlot(b1_b2_all, features = c("Ptprc", "Cd3e", "Ms4a7", "Pf4", "Plac8", "Ly6c2",
                                    "Flt3", "Itgax", "H2-Ab1", "Gata3", "Klrb1c", "Cd79a"), ncol = 4)

FeaturePlot(b1_b2_all, features = c("Apoe",
  "C1qa", "C1qb", "C1qc",
  "Pf4", "Ms4a7", "Lgmn",
  "Cd81", "Aif1", "Fcrls"), ncol = 4)

FeaturePlot(b1_b2_all, features = c("Ptprc", "Col1a1", "Pecam1", "Mbp", "Acta2", "Kcnj8", "Ncam1", "Tph1", "Mt3"))

FeaturePlot(b1_b2_all, features = "Cd3e", label = T)

###

all_immune <- subset(b1_b2_all, idents = c("3", "4", "6", "12", "14", "17", "19",
                     "20", "23", "25", "29", "31", "33", "38"))

all_immune <- NormalizeData(object = all_immune)
all_immune <- FindVariableFeatures(object = all_immune)
#length(x = VariableFeatures(object = all_immune))
all_immune <- ScaleData(object = all_immune)

###PCA
all_immune <- RunPCA(object = all_immune)
print(x = all_immune[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = all_immune, dims = 1:6)

DimPlot(object = all_immune)
all_immune <- ProjectDim(object = all_immune)

##JackStraw
all_immune <- JackStraw(object = all_immune, num.replicate = 100)
all_immune <- ScoreJackStraw(object = all_immune, dims = 1:20)
JackStrawPlot(object = all_immune, dims = 1:20)

ElbowPlot(object = all_immune, ndims = 20)

###
all_immune <- FindNeighbors(object = all_immune, dims = 1:5)
all_immune <- FindClusters(object = all_immune, resolution = 0.7)
##
all_immune <- RunUMAP(all_immune, dims = 1:5)

FeaturePlot(all_immune, features = "Ptprc", label = T, split.by = "Age")

qsave(all_immune, "/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/timecourse_all_immune.rds")

FeaturePlot(all_immune, features = c("Ptprc", "Cd3e", "Ms4a7", "Pf4", "Plac8", "Ly6c2",
                                    "Flt3", "Itgax", "H2-Ab1", "Gata3", "Klrb1c", "Cd79a"), ncol = 4)


FeaturePlot(all_immune, features = "Tcrg-C1", split.by = "Age", label = T)

Cluster5_7 <- subset(all_immune, idents = c("5", "7"))

Cluster5_7 <- NormalizeData(object = Cluster5_7)
Cluster5_7 <- FindVariableFeatures(object = Cluster5_7)
#length(x = VariableFeatures(object = all_immune))
Cluster5_7 <- ScaleData(object = Cluster5_7)

###PCA
Cluster5_7 <- RunPCA(object = Cluster5_7)
print(x = Cluster5_7[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = Cluster5_7, dims = 1:6)

DimPlot(object = Cluster5_7)
Cluster5_7 <- ProjectDim(object = Cluster5_7)

##JackStraw
Cluster5_7 <- JackStraw(object = Cluster5_7, num.replicate = 100)
Cluster5_7 <- ScoreJackStraw(object = Cluster5_7, dims = 1:20)
JackStrawPlot(object = Cluster5_7, dims = 1:20)

ElbowPlot(object = Cluster5_7, ndims = 20)

###
Cluster5_7 <- FindNeighbors(object = Cluster5_7, dims = 1:4)
Cluster5_7 <- FindClusters(object = Cluster5_7, resolution = 0.7)
##
Cluster5_7 <- RunUMAP(Cluster5_7, dims = 1:4)

FeaturePlot(Cluster5_7, features = "Tyrobp", label = T, split.by = "Age")

DimPlot(Cluster5_7, split.by = "Age", label = T)

Cluster5_7markers <- FindAllMarkers(Cluster5_7, only.pos = TRUE, min.pct = 0.25)

Cluster5_7_markers_top50 <- Cluster5_7markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.csv(Cluster5_7_markers_top50, "/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/TimeCourse_Cluster5_7_markers_top50.csv")

devtools::install_github(repo = "samuel-marsh/scCustomize")
