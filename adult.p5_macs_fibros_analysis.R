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

adult.p5_fibros <- readRDS("adult.p5_fibros.rds")
adult.p5_macs <- readRDS("adult.p5_macs.rds")

FeaturePlot_scCustom(adult.p5_fibros, features = "Spp1", label = T)

DoHeatmap(adult.p5_macs, label = T)

Idents(adult.p5_fibros) <- "Age"
adult.p5_fibros.markers <- FindMarkers(adult.p5_fibros, ident.1 = "P5", only.pos = T)



write.csv(adult.p5_fibros.markers, "adult.p5_P5positive_fibros_markers.csv")


FeaturePlot_scCustom(adult.p5_fibros, features = "Spp1", label = T)
# Eln
# Igfbp2
# Meg3
# Vcan
# Mfap2
# Dlk1
# Fn1
# Mmp13
# Spp1
# Postn
fibro_p5_markers <- c("Eln", "Igfbp2", "Meg3", "Vcan",
                      "Mfap2", "Dlk1", "Fn1", "Mmp13",
                      "Spp1", "Postn")


adult.p5_adultPositive_fibros.markers <- FindMarkers(adult.p5_fibros, ident.1 = "Adult", only.pos = T)

write.csv(adult.p5_adultPositive_fibros.markers, "adult.p5_adultPositive_fibros.markers")

FeaturePlot_scCustom(adult.p5_fibros, features = "Fxyd5", label = T)
adult_fibro_markers <- c("Gpx3", "Gpx3", "Emb", "Comp",
                         "Fgl2", "Slc47a1", "Myoc", "Ttr",
                         "Cpxm2", "Sfrp4")

# Gpx3
# Gpx3
# Comp
# Emb
# Fgl2
# Slc47a1
# Myoc
# Ttr
# Cpxm2
# Sfrp4

adult.p5_fibros_all_markers <- FindAllMarkers(object = adult.p5_fibros, only.pos = T, group.by = "Age")

adult.p5_fibros_top10_markers <- Extract_Top_Markers(marker_dataframe = adult.p5_fibros_all_markers, num_genes = 10, named_vector = FALSE,
                                    make_unique = TRUE)

Clustered_DotPlot(seurat_object = adult.p5_fibros, features = adult.p5_fibros_top10_markers, exp_color_min = -1, exp_color_max = 2,
                  colors_use_exp = PurpleAndYellow())

DoHeatmap(adult.p5_fibros, features = adult.p5_fibros_top10_markers)

DimPlot_scCustom(adult.p5_fibros, group.by = "Sex")

all_cell <- readRDS("adult_p5_all_cells.rds")


DimPlot_scCustom(all_cell, label = T, group.by = "cell_type", split.by = "Age", split_seurat = T) +
  theme(legend.text=element_text(size=15))
ggsave("allCellsp.t_adult_dimplot_celltype.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)


Idents(all_cell) <- "cell_type"
markers_perCell_byAge <- FindAllMarkers(all_cell, group.by = "cell_type", only.pos = T)


macs_fibros <- readRDS("adult.p5_macs.fibros.rds")
rm(macs_fibros)

Idents(adult.p5_macs) <- "Age"
adult.p5_adultPositive_mac.markers <- FindMarkers(adult.p5_macs, ident.1 = "Adult", only.pos = T)

adult.p5_adult_Positive_mac.markers <- write.csv(adult.p5_adultPositive_mac.markers, "adult.p5_adultPositive_mac.markers.csv")
adult.p5_adultPositive_mac.markers
FeaturePlot_scCustom(adult.p5_macs, features = "H2-Ab1", label = T)
##Mac Markers Adult
#H2-Ab1
#Cd74

adult.p5_adultPositive_mac.markers %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(adult.p5_macs, features = top10_genes_adult_mac) + NoLegend()

top10_genes_adult_mac <- c("Cd81", "Mgl2", "H2-DMa",
                     "H2-DMb1", "H2-Ab1", "H2-Aa", 
                     "H2-Eb1", "Cd74", "Ccl8")

adult.p5_P5Positive_mac.markers <- FindMarkers(adult.p5_macs, ident.1 = "P5", only.pos = T)

adult.p5_P5_Positive_mac.markers <- write.csv(adult.p5_P5Positive_mac.markers, "adult.p5_P5Positive_mac.markers.csv")
adult.p5_P5_Positive_mac.markers
FeaturePlot_scCustom(adult.p5_macs, features = "H2-Ab1", label = T)

adult.p5_P5Positive_mac.markers %>%
  top_n(n = 10, wt = avg_log2FC) -> top10_p5_macs

Idents(macs_fibros) <- "Age"
DoHeatmap(adult.p5_fibros, features = p5.adult_fibro_markers, group.by = "Age") + NoLegend()
ggsave("adult.p5_fibros_heatmap.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)


top10_genes_p5_mac <- c("Slc40a1", "Rbm3", "Hbb-bs", "Hmox1",
                           "Ednrb", "Lgals1", "Hgsnat", "Hba-a1",
                           "Hba-a2", "Hist1h2ap")

p5.adult_macs_markers <- c("Cd81", "Mgl2", "H2-DMa",
                           "H2-DMb1", "H2-Ab1", "H2-Aa", 
                           "H2-Eb1", "Cd74", "Ccl8", "Slc40a1", "Rbm3", "Hbb-bs", "Hmox1",
                           "Ednrb", "Lgals1", "Hgsnat", "Hba-a1",
                           "Hba-a2", "Hist1h2ap")

p5.adult_fibro_markers <- c("Gpx3", "Gpx3", "Emb", "Comp",
                            "Fgl2", "Slc47a1", "Myoc", "Ttr",
                            "Cpxm2", "Sfrp4", "Eln", "Igfbp2", "Meg3", "Vcan",
                            "Mfap2", "Dlk1", "Fn1", "Mmp13",
                            "Spp1", "Postn")

p5.adult_macs.fibros_markers <- c("Slc40a1", "Rbm3", "Hbb-bs", "Hmox1",
                                  "Ednrb", "Lgals1", "Hgsnat", "Hba-a1",
                                  "Hba-a2", "Hist1h2ap", "Cd81", "Mgl2", "H2-DMa",
                                  "H2-DMb1", "H2-Ab1", "H2-Aa", 
                                  "H2-Eb1", "Cd74", "Ccl8", "Gpx3", "Gpx3", "Emb", "Comp",
                                  "Fgl2", "Slc47a1", "Myoc", "Ttr",
                                  "Cpxm2", "Sfrp4", "Eln", "Igfbp2", "Meg3", "Vcan",
                                  "Mfap2", "Dlk1", "Fn1", "Mmp13",
                                  "Spp1", "Postn")

Idents(macs_fibros) <- "cell_type"
DimPlot_scCustom(macs_fibros, group.by = "Age")

FeaturePlot(macs_fibros, features = "Notch9", label = T, split.by = "Age")
