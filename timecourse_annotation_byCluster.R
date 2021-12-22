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
library(iTALK)
library(cli)



setwd("C:/Users/Stevens Lab/Desktop/LanserLIGER/data")

Sys.setenv('R_MAX_VSIZE'=memory.size(max=TRUE))

time_course_all <- readRDS("timecourse_meninges_B1_B2_filtered.rds")

time_course_all <- NormalizeData(object = time_course_all)
time_course_all <- FindVariableFeatures(object = time_course_all)
#length(x = VariableFeatures(object = all_immune))
time_course_all <- ScaleData(object = time_course_all)

###PCA
time_course_all <- RunPCA(object = time_course_all)
print(x = time_course_all[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = time_course_all, dims = 1:6)

DimPlot(object = time_course_all)
time_course_all <- ProjectDim(object = time_course_all)

##JackStraw
time_course_all <- JackStraw(object = time_course_all, num.replicate = 100)
time_course_all <- ScoreJackStraw(object = time_course_all, dims = 1:20)
JackStrawPlot(object = time_course_all, dims = 1:20)

ElbowPlot(object = time_course_all, ndims = 20)

###
time_course_all <- FindNeighbors(object = time_course_all, dims = 1:15)
time_course_all <- FindClusters(object = time_course_all, resolution = 1.2)
##
time_course_all <- RunUMAP(time_course_all, dims = 1:15)


Idents(time_course_all) = "seurat_clusters"
FeaturePlot(time_course_all, features = "Ptprc", label = T)

saveRDS(time_course_all, "timecourse_meninges_B1_B2_filtered.rds")

immune_clusters <- c("2", "4", "5", "7", "16", "18", "19", "20", "21",
                     "28")

all_immune <- subset(time_course_all, ident = c("2", "4", "5", "7", "16", "18", "19", "20", "21",
                                                   "28"))


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


ElbowPlot(object = all_immune, ndims = 20)

###
all_immune <- FindNeighbors(object = all_immune, dims = 1:7)
all_immune <- FindClusters(object = all_immune, resolution = 1.2)
##
all_immune <- RunUMAP(all_immune, dims = 1:7)

DimPlot(all_immune, label = T, split.by = "Age")

saveRDS(all_immune, "timecourse_all_immune.rds")




##Macs
Mac_score <- c("Apoe",
               "C1qa",
               "C1qb",
               "C1qc",
               "Ms4a7",
               "Lgmn",
               "Pf4",
               "Trf",
               "Hexb",
               "Ctsc")


FeaturePlot(all_immune, features = Mac_score, label = T)

all_immune <- AddModuleScore(all_immune,
                             features = list(Mac_score),
                             name="mac_enriched")

# Plot scores
FeaturePlot(all_immune, features = "mac_enriched1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

Macs <- subset(all_immune, ident = c("0", "1", "4", "5", "6", "7", "9", "11", "12",
                                     "13", "15", "16","17", "19"))



##T Cells
Tcell_score <- c("Cd3e", "Cd2", "Cd5", "Cd6",
                 "Thy1", "Cd28")
all_immune <- AddModuleScore(all_immune,
                             features = list(Tcell_score),
                             name="Tcell_enriched")

# Plot scores
FeaturePlot(all_immune,
            features = "Tcell_enriched1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

TCells <- subset(all_immune, ident = c("3"))



##Neutrophils
neutro_score <- c("S100a9", "Ly6g", "Itgam", "Chil3",
                  "Csf3r", "Camp", "Cxcr2", "Ltf")

all_immune <- AddModuleScore(all_immune,
                             features = list(neutro_score),
                             name="neutro_enriched")

# Plot scores
FeaturePlot(all_immune,
            features = "neutro_enriched1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

neutros <- subset(all_immune, ident = c("14", "20", "23"))

FeaturePlot(all_immune, features = "Il5ra")

##DCs
DCcell_score <- c("Itgax", "Lgals3", "Clec9a", "Itgae",
"Xcr1", "Itgam", "Ass1", "Batf3", "Clec10a")



all_immune <- AddModuleScore(all_immune,
                             features = list(DCcell_score),
                             name="DC_enriched")

# Plot scores
FeaturePlot(all_immune,
            features = "DC_enriched1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

DCs <- subset(all_immune, ident = c("8", "10", "15", "16", "21", "23"))


#ILC2s
ILC2_score <- c("Gata3", "Il7r")

all_immune <- AddModuleScore(all_immune,
                             features = list(ILC2_score),
                             name="ILC2_enriched")

# Plot scores
FeaturePlot(all_immune,
            features = "ILC2_enriched1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

ILC2s <- subset(all_immune, ident = c("2"))


ILC2 <- subset(all_immune, subset = ILC2_enriched1 > 0.5)


Macs$cell_type <- "Macrophages"
Tcells$cell_type <- "T Cells"
neutros$cell_type <- "Neutrophils"
DCs$cell_type <- "DCs"
ILC2$cell_type <- "ILC2s"

all_immune_labeled <- merge(x = Macs, y = c(Tcells, neutros, DCs, ILC2))


all_immune_labeled <- NormalizeData(object = all_immune_labeled)
all_immune_labeled <- FindVariableFeatures(object = all_immune_labeled)
#length(x = VariableFeatures(object = all_immune))
all_immune_labeled <- ScaleData(object = all_immune_labeled)

###PCA
all_immune_labeled <- RunPCA(object = all_immune_labeled)
print(x = all_immune_labeled[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = all_immune_labeled, dims = 1:6)

DimPlot(object = all_immune_labeled)
all_immune_labeled <- ProjectDim(object = all_immune_labeled)


ElbowPlot(object = all_immune_labeled, ndims = 20)

###
all_immune_labeled <- FindNeighbors(object = all_immune_labeled, dims = 1:12)
all_immune_labeled <- FindClusters(object = all_immune_labeled, resolution = 1.2)
##
all_immune_labeled <- RunUMAP(all_immune_labeled, dims = 1:12)



Idents(all_immune_labeled) = "cell_type"
DimPlot(all_immune_labeled, label = T, split.by = "Age")

saveRDS(all_immune_labeled, "timecourse_all_immune_labeled.rds")

all_immune_markers <- FindAllMarkers(all_immune, only.pos = T)
write.csv(all_immune_markers, "2_all_immune_markers.csv")


FeaturePlot(all_immune, features = "Cd3e")

all_macs <- subset(all_immune_labeled, ident = "Macrophages")

all_macs <- NormalizeData(object = all_macs)
all_macs <- FindVariableFeatures(object = all_macs)
#length(x = VariableFeatures(object = all_immune))
all_macs <- ScaleData(object = all_macs)

###PCA
all_macs <- RunPCA(object = all_macs)
print(x = all_macs[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = all_macs, dims = 1:6)

DimPlot(object = all_macs)
all_macs <- ProjectDim(object = all_macs)


ElbowPlot(object = all_macs, ndims = 20)

###
all_macs <- FindNeighbors(object = all_macs, dims = 1:12)
all_macs <- FindClusters(object = all_macs, resolution = 1.2)
##
all_macs <- RunUMAP(all_macs, dims = 1:12)

DimPlot(all_macs, split.by = "Age")

Idents(all_macs) = "Age"
all_mac_markers <- FindAllMarkers(all_macs, only.pos = T)
all_mac_markers_top20 <- all_mac_markers %>%
  group_by(cluster) %>% 
  top_n(20, avg_log2FC)
write.csv(all_mac_markers_top20, "all_macs_ByAge_markers.csv")


DoHeatmap(all_macs, features = all_mac_markers_top20$gene)+ scale_fill_viridis()
ggsave("Timecoursse_Heat_MoMac_top10.jpeg", device = "jpeg", 
       width = 8, height = 6, units = "in", dpi = 600)

saveRDS(all_macs, "timecourse_all_macs.rds")

all_DCs <- subset(all_immune_labeled, ident = "DCs")

all_DCs <- NormalizeData(object = all_DCs)
all_DCs <- FindVariableFeatures(object = all_DCs)
#length(x = VariableFeatures(object = all_immune))
all_DCs <- ScaleData(object = all_DCs)

###PCA
all_DCs <- RunPCA(object = all_DCs)
print(x = all_DCs[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = all_DCs, dims = 1:6)

DimPlot(object = all_DCs)
all_DCs <- ProjectDim(object = all_DCs)


ElbowPlot(object = all_DCs, ndims = 20)

###
all_DCs <- FindNeighbors(object = all_DCs, dims = 1:12)
all_DCs <- FindClusters(object = all_DCs, resolution = 1.2)
##
all_DCs <- RunUMAP(all_DCs, dims = 1:12)

DimPlot(all_DCs, split.by = "Age")

Idents(all_DCs) = "Age"
all_DCs_markers <- FindAllMarkers(all_DCs, only.pos = T)
all_DCs_markers_top20 <- all_DCs_markers %>%
  group_by(cluster) %>% 
  top_n(20, avg_log2FC)
write.csv(all_DCs_markers_top20, "all_DCs_ByAge_markers.csv")


DoHeatmap(all_DCs, features = all_DCs_markers_top20$gene)+ scale_fill_viridis()
ggsave("Timecoursse_Heat_DCs_top10.jpeg", device = "jpeg", 
       width = 8, height = 6, units = "in", dpi = 600)

saveRDS(all_DCs, "timecourse_all_DCs.rds")
