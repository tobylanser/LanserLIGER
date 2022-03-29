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

DimPlot(all_ages.cells, split.by = "Age", label = T)

Idents(all_ages.cells) <- "cell_type"
FeaturePlot_scCustom(all_ages.cells, features = "Tslp", split.by = "Age", label = T)

setwd("C:/Users/Stevens Lab/Desktop/LanserLIGER/data")

Sys.setenv('R_MAX_VSIZE'=memory.size(max=TRUE))
# install.packages("rlang")
# devtools::install_github(repo = "samuel-marsh/scCustomize")
# install.packages('NMF')
# devtools::install_github("jokergoo/circlize")
# devtools::install_github("jokergoo/ComplexHeatmap")
# devtools::install_github("sqjin/CellChat")
# install.packages("magrittr")
p5 <- readRDS("timecourse_all_P5_labeled.rds")
# 
# 
# adult <- readRDS("timecourse_all_Adult_Aged_labeled.rds")
# 
# 
# DefaultAssay(object = adult) <- "RNA"
# # adult <- SCTransform(adult, method = "glmGamPoi", vars.to.regress = "percent.mito")
# adult <- NormalizeData(adult)
# adult <- FindVariableFeatures(adult)
# adult <- ScaleData(adult)
# ###PCA
# adult <- RunPCA(object = adult, verbose = F)
# ElbowPlot(object = adult, ndims = 50)
# ###
# adult <- FindNeighbors(object = adult, dims = 1:40)
# adult <- FindClusters(object = adult, resolution = 0.8)
# adult <- RunUMAP(adult, dims = 1:40)
# 
# DimPlot_scCustom(adult, label = T, raster = F, group.by = "cell_type")
# 
# saveRDS(adult, "timecourse_adult_labeled.rds")

adult <- readRDS("timecourse_adult_labeled.rds")

all_ages.cells <- readRDS("timecourse_allAges_allCells.rds")

DimPlot_scCustom(all_ages.cells, label = T, group.by = "cell_type")
ggsave("allCells_dimplot_celltype.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)

DimPlot_scCustom(all_ages.cells, label = T, group.by = "cell_type")
ggsave("allCells_cellTypes_dimplot_celltype.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)



DimPlot_scCustom(adult, label = T, group.by = "cell_type")
ggsave("Adult_dimplot_celltype.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)

## CellChat for Adult
adult.cellchat <- createCellChat(object = adult, group.by = "cell_type")
adult.cellchatDB <- CellChatDB.mouse
adult.cellchatDB.use <- subsetDB(adult.cellchatDB, search = "Secreted Signaling")
adult.cellchat@DB <- adult.cellchatDB.use
adult.cellchat <- subsetData(adult.cellchat)

adult.cellchat <- identifyOverExpressedGenes(adult.cellchat)
adult.cellchat <- identifyOverExpressedInteractions(adult.cellchat)

adult.cellchat <- computeCommunProb(adult.cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
adult.cellchat <- filterCommunication(adult.cellchat, min.cells = 10)

adult.cellchat <- computeCommunProbPathway(adult.cellchat)
adult.cellchat <- aggregateNet(adult.cellchat)

dev.off()

groupSize <- as.numeric(table(adult.cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(adult.cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(adult.cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- adult.cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(adult.cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)


dev.off()

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(adult.cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(adult.cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(adult.cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

netAnalysis_contribution(adult.cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(adult.cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(adult.cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(adult.cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

dev.off()

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(adult.cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(adult.cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(adult.cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CCL","CXCL"),legend.pos.x = 8)
#> Note: The second link end is drawn out of sector 'CXCR4 '.
#> Note: The first link end is drawn out of sector 'CXCL12 '.

plotGeneExpression(adult.cellchat, signaling = "CXCL")

plotGeneExpression(adult.cellchat, signaling = "CXCL", enriched.only = FALSE)

# Compute the network centrality scores
adult.cellchat <- netAnalysis_computeCentrality(adult.cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(adult.cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(adult.cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(adult.cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(adult.cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(adult.cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(adult.cellchat, signaling = c("CXCL", "CCL"))

# options(PACKAGE_MAINFOLDER="C:/Users/...")
# 
# usethis::edit_r_profile()
# require(devtools)
# write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron", append = TRUE)
# packageurl <- "http://cran.us.r-project.org/src/contrib/Archive/NMF/NMF_0.22.0.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
# 
# install_version("NMF", version = "0.22.0", repos = "http://cran.us.r-project.org")
# write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron", append = TRUE)

library(NMF)
library(ggalluvial)
selectK(adult.cellchat, slot.name = "netP",
        pattern = c("outgoing", "incoming"))

nPatterns = 3
adult.cellchat <- identifyCommunicationPatterns(adult.cellchat, pattern = "outgoing", k = nPatterns)

selectK(adult.cellchat, pattern = "incoming")

usethis::edit_r_environ() 

# Set mainfolder for PACKAGE package
options(PACKAGE_MAINFOLDER="C:/Users/...")
nPatterns = 4
adult.cellchat <- identifyCommunicationPatterns(adult.cellchat, pattern = "incoming", k = nPatterns)

adult.cellchat <- computeNetSimilarity(adult.cellchat, type = "functional")
adult.cellchat <- netEmbedding(adult.cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
adult.cellchat <- netClustering(adult.cellchat, type = "functional")

adult.cellchat <- computeNetSimilarity(adult.cellchat, type = "structural")
adult.cellchat <- netEmbedding(adult.cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
adult.cellchat <- netClustering(adult.cellchat, type = "structural")

saveRDS(adult.cellchat, file = "adult.cellchat.rds")














