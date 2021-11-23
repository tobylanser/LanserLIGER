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
# remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)

setwd("C:/Users/Stevens Lab/Desktop/LanserLIGER")

Sys.setenv('R_MAX_VSIZE'=1.9e+11)

DuraCalcStrom <- readRDS("Dura_calv_stromal_2.rds")
# StromFibMyo <- readRDS("StromFibroMyo_clustered_TissueRenamed.rds")
mouseSSfibro <- readRDS("Mouse_SS_Fibro.RDS")

table(mouseSSfibro$ClustName)
mouseSSfibro$ClustName[, 1:4]


#rename project.name
DuraCalcStrom@project.name <- "DuraCalcStrom"
StromFibMyo@project.name <- "TabulaMuris"
mouseSSfibro@project.name <- "FibroblastSet"

mouseSSfibro$orig.ident <- "Turley"
StromFibMyo$orig.ident <- "Tabula Muris"
DuraCalcStrom$orig.ident <- "Meninges"

DuraCalcStrom <- UpdateSeuratObject(object = DuraCalcStrom)
mouseSSfibro <- UpdateSeuratObject(object = mouseSSfibro)
StromFibMyo <- UpdateSeuratObject(object = StromFibMyo)

fibro.big <- merge(c(DuraCalcStrom, mouseSSfibro, StromFibMyo), add.cell.ids = c("DuraCalcStrom", "Shannon Turley", "TabulaMuris"), project = "MergedFibro")

Liger3xtake2 <- seuratToLiger(list(DuraCalcStrom, StromFibMyo, mouseSSfibro), names =c("DuraCalcStrom","TabulaMuris", "Shannon Turley"))


saveRDS(fibro.big, "fibroMerged.rds")

fibro.big <- readRDS("fibroMerged.rds")

fibro.big <- NormalizeData(fibro.big)
fibro.big <- FindVariableFeatures(fibro.big)
fibro.big <- ScaleData(fibro.big, split.by = "Tissue")
#Ksuggest <- suggestK(fibro.big, num.cores = 24)
fibro.big <- RunOptimizeALS(fibro.big, k = 30, lambda = 5, split.by = "Tissue")
fibro.big <- RunQuantileNorm(fibro.big, split.by = "Tissue")
# You can optionally perform Louvain fibro.biging (`FindNeighbors` and `Findfibro.bigs`) after
# `RunQuantileNorm` according to your needs
# ElbowPlot(fibro.big, ndims =50)
seurat3x <- FindNeighbors(seurat3x, reduction = "inmf", dims = 1:ncol(Embeddings(seurat3x, "inmf"))) %>%
  FindClusters(resolution = 0.4)

# Dimensional reduction and plotting
seurat3x <- RunUMAP(seurat3x, reduction = "inmf", dims = 1:ncol(seurat3x[["inmf"]]))

FeaturePlot(seurat3x, "Pdpn", order = T)
plot1 + plot2
DimPlot(seurat3x,
        split.by ="orig.ident",
        group.by ="Tissue",
        label = TRUE,
        label.col ="white",
        label.size = 5) +DarkTheme()

DimPlot(seurat3x, raster = F, reduction = "umap", cols = c("White", "NA", "NA"), group.by = "orig.ident") +DarkTheme()
DimPlot(seurat3x, group.by = "seurat_clusters")
DimPlot(seurat3x, group.by = "Tissue")

DimPlot(seurat3x, split.by = "foxd1", group.by ="seurat_clusters", label = TRUE,
        label.size = 5, label.col ="white") +DarkTheme()

tissueplot <- DimPlot(seurat3x, label = TRUE,
        label.size = 5, group.by = "Tissue")
clusterplot <- DimPlot(seurat3x, label = TRUE,
        label.size = 5, group.by = "seurat_clusters")
# DimPlot(seurat3x, label = TRUE,
#         label.size = 5, group.by = "Tissue", split.by = "orig.ident")
clusterplot
tissueplot + clusterplot
features <- c("Bmp4", "Ccl19", "Coch", "Cxcl12", "Col15a1", "Comp", "Fbln1", "Hhip", "Npnt", "Pi16", "Dpt")

FeaturePlot(DuraCalcStrom, label = TRUE,
            label.size = 5, features = features)

FeaturePlot(DuraCalcStrom, label = TRUE,
            label.size = 5, features = "Dpt")

DoHeatmap(DuraCalcStrom, group.by = "Tissue", features = features) + DarkTheme()

DimPlot(DuraCalcStrom, group.by = "Tissue") + DarkTheme()

library(RColorBrewer)
DarkThemeFeature + DarkTheme()
FeaturePlot(seurat3x, features = c("Foxd1"), raster = F, order = T, cols = c("black","white")) + DarkTheme()

plot + theme(panel.background = element_rect(fill = 'green', colour = 'red'))

VlnPlot(seurat3x, group.by ="Tissue", features ="Gm42418") + DarkTheme()
DotPlot(seurat3x, group.by = "Tissue", features = "Foxd1", cols = c("black","white")) + RotatedAxis() + DarkTheme()


pathway_foxd1_plot <- DEenrichRPlot(seurat3x, ident.1 = foxd1, max.genes = 1000, enrich.database = "WikiPathways_2019_Mouse")


# Get number of cells per cluster and per sample of origin
table(seurat3x$RNA_snn_res.0.2, seurat3x$orig.ident)

head(AverageExpression(seurat3x, features = "Foxd1", group.by = "Tissue"))

rm(rna_assays)
#pathway analysis
install.packages('enrichR')
Idents(object = seurat3x) <- seurat3x@meta.data$"orig.ident"
DEPlot <- DEenrichRPlot(seurat3x, ident.1 = "orig.ident", ident.2 = "Tissue", enrich.database ="WikiPathways_2019_Mouse", max.genes = 1000)


devtools::install_github(repo = "samuel-marsh/scCustomize", ref = "develop", auth_token = "03a4982b27ad3b4304c6d5c3fb5da45a9106368d")
library(scCustomize)





tibble(
  cluster = seurat3x$seurat_clusters,
  Dataset = seurat3x$foxd1
) %>%
  group_by(cluster, Dataset) %>%
  count() %>%
  group_by(cluster) %>%
  mutate(
    percent=(100*n)/sum(n)
  ) %>%
  ungroup() %>%
  mutate(
    cluster=paste("Cluster", cluster)
  ) %>%
  ggplot(aes(x="",y=percent, fill=Dataset)) +
  geom_col(width=1) +
  coord_polar("y", start=0) +
  facet_wrap(vars(cluster)) +  
  theme(axis.text.x=element_blank()) +
  xlab(NULL) +
  ylab(NULL)



lapply(
  levels(fibro.big[["seurat_clusters"]][[1]]),
  function(x)FindMarkers(fibro.big,ident.1 = x,min.pct = 0.25)
) -> cluster.markers


# This simply adds the cluster number to the results of FindMarkers
sapply(0:(length(cluster.markers)-1),function(x) {
  cluster.markers[[x+1]]$gene <<- rownames(cluster.markers[[x+1]])
  cluster.markers[[x+1]]$cluster <<- x
})




cluster.markers %>%
  group_by(cluster) %>%
  slice(1) %>%
  pull(gene) -> best.wilcox.gene.per.cluster


fibro.big2 <- readRDS("fibroMerged.rds")
# split the dataset into a list of two seurat objects (stim and CTRL)
fibro.big2.list <- SplitObject(fibro.big2, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
fibro.big2.list <- lapply(X = fibro.big2.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = fibro.big2.list)

fibro.anchors <- FindIntegrationAnchors(object.list = fibro.big2.list, anchor.features = features)

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = fibro.anchors)




