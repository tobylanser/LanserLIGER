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

ILC_Tcells <- qread("/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/rds_files/ILC_Tcells_timecourse.rds")


###Aged ILC2s & T Cells 
Idents(ILC_Tcells) <- "Age"
aged_ILC_Tcells <- subset(ILC_Tcells, idents = "Aged")

aged_ILC_Tcells <- NormalizeData(object = aged_ILC_Tcells)
aged_ILC_Tcells <- FindVariableFeatures(object = aged_ILC_Tcells)
aged_ILC_Tcells <- ScaleData(object = aged_ILC_Tcells)

###PCA
aged_ILC_Tcells <- RunPCA(object = aged_ILC_Tcells)
print(x = aged_ILC_Tcells[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = aged_ILC_Tcells, dims = 1:6)

DimPlot(object = aged_ILC_Tcells)
aged_ILC_Tcells <- ProjectDim(object = aged_ILC_Tcells)

##JackStraw
aged_ILC_Tcells <- JackStraw(object = aged_ILC_Tcells, num.replicate = 100)
aged_ILC_Tcells <- ScoreJackStraw(object = aged_ILC_Tcells, dims = 1:20)
JackStrawPlot(object = aged_ILC_Tcells, dims = 1:20)

ElbowPlot(object = aged_ILC_Tcells, ndims = 20)

###
aged_ILC_Tcells <- FindNeighbors(object = aged_ILC_Tcells, dims = 1:5)
aged_ILC_Tcells <- FindClusters(object = aged_ILC_Tcells, resolution = 0.7)
##
aged_ILC_Tcells <- RunUMAP(aged_ILC_Tcells, dims = 1:5)

FeaturePlot(aged_ILC_Tcells, features = "Cd3e", label = T, order = T)

DimPlot_scCustom(aged_ILC_Tcells, label = T, split_seurat = T) 

qsave(aged_ILC_Tcells, "/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/rds_files/Aged_ILC_Tcells_timecourse.rds")

### Aged Fibroblasts

b1_b2_all_nonimmune <- qread("timecourse_all_nonimmune.rds")

Idents(b1_b2_all_nonimmune) <- "Age"
aged_b1_b2_all_nonimmune <- subset(b1_b2_all_nonimmune, idents = "Aged")

aged_b1_b2_all_nonimmune <- NormalizeData(object = aged_b1_b2_all_nonimmune)
aged_b1_b2_all_nonimmune <- FindVariableFeatures(object = aged_b1_b2_all_nonimmune)
aged_b1_b2_all_nonimmune <- ScaleData(object = aged_b1_b2_all_nonimmune)

###PCA
aged_b1_b2_all_nonimmune <- RunPCA(object = aged_b1_b2_all_nonimmune)
print(x = aged_b1_b2_all_nonimmune[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = aged_b1_b2_all_nonimmune, dims = 1:6)

DimPlot(object = aged_b1_b2_all_nonimmune)
aged_b1_b2_all_nonimmune <- ProjectDim(object = aged_b1_b2_all_nonimmune)

##JackStraw
aged_b1_b2_all_nonimmune <- JackStraw(object = aged_b1_b2_all_nonimmune, num.replicate = 100)
aged_b1_b2_all_nonimmune <- ScoreJackStraw(object = aged_b1_b2_all_nonimmune, dims = 1:20)
JackStrawPlot(object = aged_b1_b2_all_nonimmune, dims = 1:20)

ElbowPlot(object = aged_b1_b2_all_nonimmune, ndims = 20)

###
aged_b1_b2_all_nonimmune <- FindNeighbors(object = aged_b1_b2_all_nonimmune, dims = 1:8)
aged_b1_b2_all_nonimmune <- FindClusters(object = aged_b1_b2_all_nonimmune, resolution = 0.7)
##
aged_b1_b2_all_nonimmune <- RunUMAP(aged_b1_b2_all_nonimmune, dims = 1:8)

FeaturePlot(aged_b1_b2_all_nonimmune, features = "Loxl1", label = T, order = T)

Idents(aged_b1_b2_all_nonimmune) <- "seurat_clusters"
aged_fibros <- subset(aged_b1_b2_all_nonimmune, idents = c("1", "2", "3", "11", "14"))

aged_fibros <- NormalizeData(object = aged_fibros)
aged_fibros <- FindVariableFeatures(object = aged_fibros)
aged_fibros <- ScaleData(object = aged_fibros)

###PCA
aged_fibros <- RunPCA(object = aged_fibros)
print(x = aged_fibros[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = aged_fibros, dims = 1:6)

DimPlot(object = aged_fibros)
aged_fibros <- ProjectDim(object = aged_fibros)

##JackStraw
aged_fibros <- JackStraw(object = aged_fibros, num.replicate = 100)
aged_fibros <- ScoreJackStraw(object = aged_fibros, dims = 1:20)
JackStrawPlot(object = aged_fibros, dims = 1:20)

ElbowPlot(object = aged_fibros, ndims = 20)

###
aged_fibros <- FindNeighbors(object = aged_fibros, dims = 1:5)
aged_fibros <- FindClusters(object = aged_fibros, resolution = 0.7)
##
aged_fibros <- RunUMAP(aged_fibros, dims = 1:5)

FeaturePlot(aged_fibros, features = "Col1a1", label = T, order = T)

aged_fibros$cell_type <- "Fibros"
qsave(aged_fibros, "/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/rds_files/Aged_fibros_timecourse.rds")

###iTALK on aged ILCs, T cells, and Fibroblasts

aged_ILC_Tcells_fibros <- merge(x = aged_ILC_Tcells, y = aged_fibros)

aged_ILC_Tcells_fibros <- NormalizeData(object = aged_ILC_Tcells_fibros)
aged_ILC_Tcells_fibros <- FindVariableFeatures(object = aged_ILC_Tcells_fibros)
#length(x = VariableFeatures(object = all_immune))
aged_ILC_Tcells_fibros <- ScaleData(object = aged_ILC_Tcells_fibros)

###PCA
aged_ILC_Tcells_fibros <- RunPCA(object = aged_ILC_Tcells_fibros)
print(x = aged_ILC_Tcells_fibros[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = aged_ILC_Tcells_fibros, dims = 1:6)

DimPlot(object = aged_ILC_Tcells_fibros)
aged_ILC_Tcells_fibros <- ProjectDim(object = aged_ILC_Tcells_fibros)

##JackStraw
aged_ILC_Tcells_fibros <- JackStraw(object = aged_ILC_Tcells_fibros, num.replicate = 100)
aged_ILC_Tcells_fibros <- ScoreJackStraw(object = aged_ILC_Tcells_fibros, dims = 1:20)
JackStrawPlot(object = aged_ILC_Tcells_fibros, dims = 1:20)

ElbowPlot(object = aged_ILC_Tcells_fibros, ndims = 20)

###
aged_ILC_Tcells_fibros <- FindNeighbors(object = aged_ILC_Tcells_fibros, dims = 1:4)
aged_ILC_Tcells_fibros <- FindClusters(object = aged_ILC_Tcells_fibros, resolution = 0.7)
##
aged_ILC_Tcells_fibros <- RunUMAP(aged_ILC_Tcells_fibros, dims = 1:4)

Idents(aged_ILC_Tcells_fibros) <- "cell_type"

DimPlot_scCustom(aged_ILC_Tcells_fibros, label = T)

FeaturePlot_scCustom(aged_ILC_Tcells_fibros, features = "Gata3", label = T)

df <- t(as.matrix(GetAssayData(object = aged_ILC_Tcells_fibros, slot = "counts")))

colnames(df) <- toupper(colnames(df))


ct <- as.matrix(aged_ILC_Tcells_fibros$cell_type)

colnames(ct) <- "cell_type"

merged <- cbind(df, ct)
write.table(merged, '/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/Aged_2geneXcell_Matrix_ILC_Tcells_fibros.txt', 
            sep = '\t', row.names = T, col.names = T, quote = F)



Aged_ILC_TCell_fibros_iTALK <-read.table('/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/2geneXcell_Matrix_ILC_Tcells_fibros_TC.txt', sep='\t', header=T, stringsAsFactors = F)
## highly expressed ligand-receptor pairs

# find top 50 percent highly expressed genes
highly_exprs_genes<-rawParse(Aged_ILC_TCell_fibros_iTALK,top_genes=50,stats='mean')
# find the ligand-receptor pairs from highly expressed genes
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_col<-structure(c('#4a84ad','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b','#9ef49a'),names=unique(Aged_ILC_TCell_fibros_iTALK$cell_type))
par(mfrow=c(1,2))
res<-NULL
for(comm_type in comm_list){
  res_cat<-FindLR(highly_exprs_genes,datatype='mean count',comm_type=comm_type)
  res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
  #plot by ligand category
  #overall network plot
  NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
  #top 20 ligand-receptor pairs
  LRPlot(res_cat[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_mean_exprs[1:20],link.arr.width=res_cat$cell_to_mean_exprs[1:20])
  title(comm_type)
  res<-rbind(res,res_cat)
}
res<-res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:20,]
NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
LRPlot(res[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res$cell_from_mean_exprs[1:20],link.arr.width=res$cell_to_mean_exprs[1:20])

## significant ligand-receptor pairs between compare groups

# randomly assign the compare group to each sample
Aged_ILC_TCell_fibros_iTALK<-Aged_ILC_TCell_fibros_iTALK %>% mutate(compare_group=sample(2,nrow(Aged_ILC_TCell_fibros_iTALK),replace=TRUE))
# find DEGenes of gdT cells and NK cells between these 2 groups
deg_t<-DEG(Aged_ILC_TCell_fibros_iTALK %>% filter(cell_type=='gdTcells'),method='Wilcox',contrast=c(2,1))
deg_ilc2s<-DEG(Aged_ILC_TCell_fibros_iTALK %>% filter(cell_type=='Fibros'),method='Wilcox',contrast=c(2,1))
# find significant ligand-receptor pairs and do the plotting
par(mfrow=c(1,2))
res<-NULL
for(comm_type in comm_list){
  res_cat<-FindLR(deg_t,deg_nk,datatype='DEG',comm_type=comm_type)
  res_cat<-res_cat[order(res_cat$cell_from_logFC*res_cat$cell_to_logFC,decreasing=T),]
  #plot by ligand category
  if(nrow(res_cat)==0){
    next
  }else if(nrow(res_cat>=20)){
    LRPlot(res_cat[1:20,],datatype='DEG',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_logFC[1:20],link.arr.width=res_cat$cell_to_logFC[1:20])
  }else{
    LRPlot(res_cat,datatype='DEG',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_logFC,link.arr.width=res_cat$cell_to_logFC)
  }
  NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
  title(comm_type)
  res<-rbind(res,res_cat)
}
if(is.null(res)){
  print('No significant pairs found')
}else if(nrow(res)>=20){
  res<-res[order(res$cell_from_logFC*res$cell_to_logFC,decreasing=T),][1:20,]
  NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
  LRPlot(res[1:20,],datatype='DEG',cell_col=cell_col,link.arr.lwd=res$cell_from_logFC[1:20],link.arr.width=res$cell_to_logFC[1:20])
}else{
  NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
  LRPlot(res,datatype='DEG',cell_col=cell_col,link.arr.lwd=res$cell_from_logFC,link.arr.width=res$cell_to_logFC)
}
# I just randomly assigned the compare group to samples which has no biological difference for showing how to use the package.
# So there should be no significant genes to be expected.
