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
library(iTALK)
library(cli)



setwd("C:/Users/Stevens Lab/Desktop/LanserLIGER/data")

Sys.setenv('R_MAX_VSIZE'=memory.size(max=TRUE))

time_course_all <- readRDS("timecourse_meninges_B1_B2_filtered.rds")


DimPlot(time_course_all, label = T)

Idents(time_course_all) = "seurat_clusters"
FeaturePlot(time_course_all, features = c("Col1a1", "Ptprc"))
FeaturePlot(time_course_all, features = c("Ms4a1", "Vpreb1", "Il7r", "Rag2"), label = T)                                    
FeaturePlot(time_course_all, features = c("Plac8", "Ly6c2", "Flt3", "Itgax"), label = T)
FeaturePlot(time_course_all, features = c("H2-Ab1", "Gata3", "Klrb1c", "Cd79a"), label = T)



immune_score <-c("Ptprc")

time_course_all <- AddModuleScore(time_course_all,
                                features = list(immune_score),
                                name="immune_enriched")

# Plot scores
FeaturePlot(time_course_all,
            features = "immune_enriched1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

timecourse_all_immune <- subset(time_course_all, subset = immune_enriched1 > 0.01)


all_immune <- timecourse_all_immune
rm(timecourse_all_immune)
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
all_immune <- FindNeighbors(object = all_immune, dims = 1:7)
all_immune <- FindClusters(object = all_immune, resolution = 0.7)
##
all_immune <- RunUMAP(all_immune, dims = 1:7)

DimPlot(all_immune, label = T, split.by = "Age")
FeaturePlot(all_immune, features = "Ptprc")
saveRDS(all_immune, "timecourse_all_immune.rds")





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
mac_markerPlot <- FeaturePlot(all_immune,
            features = "mac_enriched1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

mac_markerPlot

Macs <- subset(all_immune, subset = mac_enriched1 > 0.5)
Macs$cell_type <- "Macrophages"


Bcell_score <- c("Cd79a", "Cd19", "Bach2", "Blnk",
                 "Ms4a1", "Vpreb1", "Il7r", "Rag2")

all_immune <- AddModuleScore(all_immune,
                                features = list(Bcell_score),
                                name="Bcell_enriched")

# Plot scores
FeaturePlot(all_immune,
            features = "Bcell_enriched1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



Bcells <- subset(all_immune, subset = Bcell_enriched1 > 0.5)
Bcells$cell_type <- "B Cells"

FeaturePlot(all_immune, features = "Cd28", label = T)
FeaturePlot(all_immune, features = c("Cd3e", "Cd2"), label = T)

Tcell_score <- c("Cd3e", "Cd2", "Cd5", "Cd6",
                 "Thy1", "Cd28")
all_immune <- AddModuleScore(all_immune,
                             features = list(Tcell_score),
                             name="Tcell_enriched")

# Plot scores
FeaturePlot(all_immune,
            features = "Tcell_enriched1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



Tcells <- subset(all_immune, subset = Tcell_enriched1 > 0.5)
Tcells$cell_type <- "T Cells"


tcells$cell_type <- "T Cells"

FeaturePlot(all_immune, features = "Rorc", label = T)
FeaturePlot(all_immune, features = c("Cd3e", "Cd2"), label = T)


neutro_score <- c("S100a9", "Ly6g", "Itgam", "Chil3",
                 "Csf3r", "Camp", "Cxcr2", "Ltf")

all_immune <- AddModuleScore(all_immune,
                             features = list(neutro_score),
                             name="neutro_enriched")

# Plot scores
FeaturePlot(all_immune,
            features = "neutro_enriched1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



neutros <- subset(all_immune, subset = neutro_enriched1 > 0.5)


neutros$cell_type <- "Neutrophils"


FeaturePlot(all_immune, features = "Sirpa", label = T, split.by = "Age")
FeaturePlot(all_immune, features = c("Cd3e", "Cd2"), label = T)


DCcell_score <- c("Itgax", "Lgals3", "Clec9a", "Itgae",
                  "Xcr1", "Itgam", "Ass1", "Batf3", "Clec10a")

FeaturePlot(all_immune, features = DCcell_score)

all_immune <- AddModuleScore(all_immune,
                             features = list(DCcell_score),
                             name="DC_enriched")

# Plot scores
FeaturePlot(all_immune,
            features = "DC_enriched1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

DC_MarkerPlot

DCs <- subset(all_immune, subset = DC_enriched1 > 0.1)

DimPlot(all_immune, label = T)

DCs$cell_type <- "DCs"


FeaturePlot(all_immune, features = "Sparc", label = T)
FeaturePlot(all_immune, features = c("Col1a1", "Ptprc"), label = T)


ILC2_score <- c("Gata3", "Il7r")

all_immune <- AddModuleScore(all_immune,
                             features = list(ILC2_score),
                             name="ILC2_enriched")

# Plot scores
FeaturePlot(all_immune,
            features = "ILC2_enriched1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



ILC2 <- subset(all_immune, subset = ILC2_enriched1 > 0.5)


ILC2$cell_type <- "ILC2s"


all_immune_labeled <- merge(x = Macs, y = c(Tcells, neutros, Bcells, DCs, ILC2))


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

##JackStraw
all_immune_labeled <- JackStraw(object = all_immune_labeled, num.replicate = 100)
all_immune_labeled <- ScoreJackStraw(object = all_immune_labeled, dims = 1:20)
JackStrawPlot(object = all_immune_labeled, dims = 1:20)

ElbowPlot(object = all_immune_labeled, ndims = 20)

###
all_immune_labeled <- FindNeighbors(object = all_immune_labeled, dims = 1:12)
all_immune_labeled <- FindClusters(object = all_immune_labeled, resolution = 1.2)
##
all_immune_labeled <- RunUMAP(all_immune_labeled, dims = 1:12)



Idents(all_immune_labeled) = "cell_type"
DimPlot(all_immune_labeled, label = T)

all_immune_markers <- FindAllMarkers(all_immune_labeled, only.pos = T)
write.csv(all_immune_markers, "all_immune_markers.csv")


                               