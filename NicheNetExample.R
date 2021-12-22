

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
library(Cairo)
library(iTALK)
library(nichenetr)
library(garnett)
library(Hmisc)
library(survival)
library(lattice)


setwd("/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/rds_files")



ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network



adult_ILC_Tcells_fibros <- qread("Adult_ILC_Tcells_fibros_timecourse.rds")


# indicated cell types should be cell class identities
# check via: 
# seuratObj %>% Idents() %>% table()
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = adult_ILC_Tcells_fibros, 
  receiver = "CD8 T", 
  condition_colname = "cell_type", condition_oi = "cell_type", condition_reference = "SS", 
  sender = c("CD4","CD8", "Fibros", "gdTcells", "ILC2s", "NKCells"), 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")
## [1] "Read in and process NicheNet's networks"
table(adult_ILC_Tcells_fibros$cell_type)
