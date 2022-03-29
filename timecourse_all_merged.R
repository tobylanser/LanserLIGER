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
library(celltalker)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)



setwd("C:/Users/Stevens Lab/Desktop/LanserLIGER/data")

Sys.setenv('R_MAX_VSIZE'=memory.size(max=TRUE))

all_immune <- readRDS("timecourse_all_immune_labeled.rds")
all_non_immune <- readRDS("timecourse_nonimmune_labeled.rds")

timecourse_all <- merge(all_immune, all_non_immune)

timecourse_all <- SCTransform(timecourse_all)
timecourse_all <- RunPCA(object = timecourse_all)

ElbowPlot(object = timecourse_all, ndims = 50)

###
timecourse_all <- FindNeighbors(object = timecourse_all, dims = 1:40)
timecourse_all <- FindClusters(object = timecourse_all, resolution = 1.2)
##
timecourse_all <- RunUMAP(timecourse_all, dims = 1:40)

Idents(timecourse_all) = "cell_type"
DimPlot_scCustom(timecourse_all, label = T, group.by = "Sex")
Aggsave("Timecourse_all_labeled.jpeg", device = "jpeg", 
       width = 10, height = 6, units = "in", dpi = 600)

saveRDS(timecourse_all, "timecourse_all_labeled.rds")
timecourse_all <- readRDS("timecourse_all_labeled.rds")

##CellChat
cellchat <- createCellChat(object = timecourse_all, group.by = "cell_type")
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
library(plan)
install.packages("parallelly")
availableCores()
# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
num_interactions <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
weighted_ints <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

