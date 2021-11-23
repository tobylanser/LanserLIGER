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

write.csv(seurat3x@reductions$umap@cell.embeddings, file = paste("UMAP_INMFCoordindates.csv"))

seurat3x@reductions$umap@cell.embeddings
umap_coord <- seurat3x@reductions$umap@cell.embeddings
write.table(umap_coord, file = paste("UMAP_INMFCoordindates.txt"), quote = F, col.names = F, row.names = T, sep = "\t")

seurat3x_active.ident <- as.data.frame(seurat3x@active.ident)

fwrite(seurat3x_active.ident, file = paste("seurat3x_active.ident.txt"), col.names = F, row.names = F, sep = "\t")
rm(rna_assays)
rna_assays <- seurat3x@assays$RNA@scale.data

#write to csv
library(data.table)
seurat3x_data_to_write_out <- as.data.frame(as.matrix(rna_assays))

fwrite(seurat3x_data_to_write_out, file = paste("seurat3xScaled_GeneXCellData.csv"), col.names = T, quote = F, row.names = T)








