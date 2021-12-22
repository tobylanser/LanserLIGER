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


DimPlot(aged_fibros, split.by = "Age")
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

qsave(aged_ILC_Tcells_fibros, "/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/rds_files/Aged_ILC_Tcells_fibros_timecourse.rds")

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

###Adult ILC2s & T Cells 
Idents(ILC_Tcells) <- "Age"
Adult_ILC_Tcells <- subset(ILC_Tcells, idents = "Adult")

Adult_ILC_Tcells <- NormalizeData(object = Adult_ILC_Tcells)
Adult_ILC_Tcells <- FindVariableFeatures(object = Adult_ILC_Tcells)
Adult_ILC_Tcells <- ScaleData(object = Adult_ILC_Tcells)

###PCA
Adult_ILC_Tcells <- RunPCA(object = Adult_ILC_Tcells)
print(x = Adult_ILC_Tcells[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = Adult_ILC_Tcells, dims = 1:6)

DimPlot(object = Adult_ILC_Tcells)
Adult_ILC_Tcells <- ProjectDim(object = Adult_ILC_Tcells)

##JackStraw
Adult_ILC_Tcells <- JackStraw(object = Adult_ILC_Tcells, num.replicate = 100)
Adult_ILC_Tcells <- ScoreJackStraw(object = Adult_ILC_Tcells, dims = 1:20)
JackStrawPlot(object = Adult_ILC_Tcells, dims = 1:20)

ElbowPlot(object = Adult_ILC_Tcells, ndims = 20)

###
Adult_ILC_Tcells <- FindNeighbors(object = Adult_ILC_Tcells, dims = 1:5)
Adult_ILC_Tcells <- FindClusters(object = Adult_ILC_Tcells, resolution = 0.7)
##
Adult_ILC_Tcells <- RunUMAP(Adult_ILC_Tcells, dims = 1:5)

FeaturePlot(Adult_ILC_Tcells, features = "Cd3e", label = T, order = T)

DimPlot_scCustom(Adult_ILC_Tcells, label = T, split_seurat = T) 

qsave(Adult_ILC_Tcells, "/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/rds_files/Adult_ILC_Tcells_timecourse.rds")

### Aged Fibroblasts

b1_b2_all_nonimmune <- qread("timecourse_all_nonimmune.rds")

DimPlot(b1_b2_all_nonimmune, split.by = "Age")
DimPlot(b1_b2_all_nonimmune)
FeaturePlot_scCustom(b1_b2_all_nonimmune, features = "Pecam1", split.by = "Age")
FeaturePlot_scCustom(b1_b2_all_nonimmune, features = "Pecam1")

Idents(b1_b2_all_nonimmune) <- "Age"
Adult_b1_b2_all_nonimmune <- subset(b1_b2_all_nonimmune, idents = "Adult")

Adult_b1_b2_all_nonimmune <- NormalizeData(object = Adult_b1_b2_all_nonimmune)
Adult_b1_b2_all_nonimmune <- FindVariableFeatures(object = Adult_b1_b2_all_nonimmune)
Adult_b1_b2_all_nonimmune <- ScaleData(object = Adult_b1_b2_all_nonimmune)

###PCA
Adult_b1_b2_all_nonimmune <- RunPCA(object = Adult_b1_b2_all_nonimmune)
print(x = Adult_b1_b2_all_nonimmune[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = Adult_b1_b2_all_nonimmune, dims = 1:6)

DimPlot(object = Adult_b1_b2_all_nonimmune)
Adult_b1_b2_all_nonimmune <- ProjectDim(object = Adult_b1_b2_all_nonimmune)

##JackStraw
Adult_b1_b2_all_nonimmune <- JackStraw(object = Adult_b1_b2_all_nonimmune, num.replicate = 100)
Adult_b1_b2_all_nonimmune <- ScoreJackStraw(object = Adult_b1_b2_all_nonimmune, dims = 1:20)
JackStrawPlot(object = Adult_b1_b2_all_nonimmune, dims = 1:20)

ElbowPlot(object = Adult_b1_b2_all_nonimmune, ndims = 20)

###
Adult_b1_b2_all_nonimmune <- FindNeighbors(object = Adult_b1_b2_all_nonimmune, dims = 1:8)
Adult_b1_b2_all_nonimmune <- FindClusters(object = Adult_b1_b2_all_nonimmune, resolution = 0.7)
##
Adult_b1_b2_all_nonimmune <- RunUMAP(Adult_b1_b2_all_nonimmune, dims = 1:8)

FeaturePlot(Adult_b1_b2_all_nonimmune, features = "Col1a1", label = T, order = T)

Idents(Adult_b1_b2_all_nonimmune) <- "seurat_clusters"
adult_fibros <- subset(Adult_b1_b2_all_nonimmune, idents = c("0", "1", "3", "5", "6", "7", "13"))

adult_fibros <- NormalizeData(object = adult_fibros)
adult_fibros <- FindVariableFeatures(object = adult_fibros)
adult_fibros <- ScaleData(object = adult_fibros)

###PCA
adult_fibros <- RunPCA(object = adult_fibros)
print(x = adult_fibros[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = adult_fibros, dims = 1:6)

DimPlot(object = adult_fibros)
adult_fibros <- ProjectDim(object = adult_fibros)

##JackStraw
adult_fibros <- JackStraw(object = adult_fibros, num.replicate = 100)
adult_fibros <- ScoreJackStraw(object = adult_fibros, dims = 1:20)
JackStrawPlot(object = adult_fibros, dims = 1:20)

ElbowPlot(object = adult_fibros, ndims = 20)

###
adult_fibros <- FindNeighbors(object = adult_fibros, dims = 1:7)
adult_fibros <- FindClusters(object = adult_fibros, resolution = 0.7)
##
adult_fibros <- RunUMAP(adult_fibros, dims = 1:7)

FeaturePlot(adult_fibros, features = "Col1a1", label = T, order = T)

adult_fibros$cell_type <- "Fibros"
qsave(adult_fibros, "/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/rds_files/Adult_fibros_timecourse.rds")

###iTALK on adult ILCs, T cells, and Fibroblasts

adultILC_Tcells_fibros <- merge(x = Adult_ILC_Tcells, y = adult_fibros)

adultILC_Tcells_fibros <- NormalizeData(object = adultILC_Tcells_fibros)
adultILC_Tcells_fibros <- FindVariableFeatures(object = adultILC_Tcells_fibros)
#length(x = VariableFeatures(object = all_immune))
adultILC_Tcells_fibros <- ScaleData(object = adultILC_Tcells_fibros)

###PCA
adultILC_Tcells_fibros <- RunPCA(object = adultILC_Tcells_fibros)
print(x = adultILC_Tcells_fibros[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = adultILC_Tcells_fibros, dims = 1:6)

DimPlot(object = adultILC_Tcells_fibros)
adultILC_Tcells_fibros <- ProjectDim(object = adultILC_Tcells_fibros)

##JackStraw
adultILC_Tcells_fibros <- JackStraw(object = adultILC_Tcells_fibros, num.replicate = 100)
adultILC_Tcells_fibros <- ScoreJackStraw(object = adultILC_Tcells_fibros, dims = 1:20)
JackStrawPlot(object = adultILC_Tcells_fibros, dims = 1:20)

ElbowPlot(object = adultILC_Tcells_fibros, ndims = 20)

###
adultILC_Tcells_fibros <- FindNeighbors(object = adultILC_Tcells_fibros, dims = 1:5)
adultILC_Tcells_fibros <- FindClusters(object = adultILC_Tcells_fibros, resolution = 0.7)
##
adultILC_Tcells_fibros <- RunUMAP(adultILC_Tcells_fibros, dims = 1:5)

qsave(adultILC_Tcells_fibros, "/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/rds_files/Adult_ILC_Tcells_fibros_timecourse.rds")

Idents(adultILC_Tcells_fibros) <- "cell_type"

DimPlot_scCustom(adultILC_Tcells_fibros, label = T)

FeaturePlot_scCustom(adultILC_Tcells_fibros, features = "Gata3", label = T)

df <- t(as.matrix(GetAssayData(object = adultILC_Tcells_fibros, slot = "counts")))

colnames(df) <- toupper(colnames(df))

all_ILC_Tcells <- qread("ILC_Tcells_timecourse.rds")

ct <- as.matrix(adultILC_Tcells_fibros$cell_type)

colnames(ct) <- "cell_type"

merged <- cbind(df, ct)
write.table(merged, '/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/Adult_2geneXcell_Matrix_ILC_Tcells_fibros.txt', 
            sep = '\t', row.names = T, col.names = T, quote = F)



Adult_ILC_TCell_fibros_iTALK <-read.table('/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/2geneXcell_Matrix_ILC_Tcells_fibros_TC.txt', sep='\t', header=T, stringsAsFactors = F)
## highly expressed ligand-receptor pairs

# find top 50 percent highly expressed genes
adult_highly_exprs_genes<-rawParse(Adult_ILC_TCell_fibros_iTALK,top_genes=50,stats='mean')
# find the ligand-receptor pairs from highly expressed genes
comm_list<-c('growth factor')
cell_col<-structure(c('#4a84ad','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b','#9ef49a'),names=unique(Adult_ILC_TCell_fibros_iTALK$cell_type))
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
Adult_ILC_TCell_fibros_iTALK<-Adult_ILC_TCell_fibros_iTALK %>% mutate(compare_group=sample(2,nrow(Adult_ILC_TCell_fibros_iTALK),replace=TRUE))
# find DEGenes of gdT cells and NK cells between these 2 groups
deg_t<-DEG(Adult_ILC_TCell_fibros_iTALK %>% filter(cell_type=='gdTcells'),method='Wilcox',contrast=c(2,1))
deg_fibros<-DEG(Adult_ILC_TCell_fibros_iTALK %>% filter(cell_type=='Fibros'),method='Wilcox',contrast=c(2,1))
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

###P5 ILC2s & T Cells 
Idents(ILC_Tcells) <- "Age"
p5_ILC_Tcells <- subset(ILC_Tcells, idents = "P5")

p5_ILC_Tcells <- NormalizeData(object = p5_ILC_Tcells)
p5_ILC_Tcells <- FindVariableFeatures(object = p5_ILC_Tcells)
p5_ILC_Tcells <- ScaleData(object = p5_ILC_Tcells)

###PCA
p5_ILC_Tcells <- RunPCA(object = p5_ILC_Tcells)
print(x = p5_ILC_Tcells[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = p5_ILC_Tcells, dims = 1:6)

DimPlot(object = p5_ILC_Tcells)
p5_ILC_Tcells <- ProjectDim(object = p5_ILC_Tcells)

##JackStraw
p5_ILC_Tcells <- JackStraw(object = p5_ILC_Tcells, num.replicate = 100)
p5_ILC_Tcells <- ScoreJackStraw(object = p5_ILC_Tcells, dims = 1:20)
JackStrawPlot(object = p5_ILC_Tcells, dims = 1:20)

ElbowPlot(object = p5_ILC_Tcells, ndims = 20)

###
p5_ILC_Tcells <- FindNeighbors(object = p5_ILC_Tcells, dims = 1:5)
p5_ILC_Tcells <- FindClusters(object = p5_ILC_Tcells, resolution = 0.7)
##
p5_ILC_Tcells <- RunUMAP(p5_ILC_Tcells, dims = 1:5)

FeaturePlot(p5_ILC_Tcells, features = "Cd3e", label = T, order = T)

DimPlot_scCustom(p5_ILC_Tcells, label = T, split_seurat = T) 

qsave(p5_ILC_Tcells, "/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/rds_files/p5_ILC_Tcells_timecourse.rds")

### Aged Fibroblasts

b1_b2_all_nonimmune <- qread("timecourse_all_nonimmune.rds")

DimPlot(b1_b2_all_nonimmune, split.by = "Age")
DimPlot(b1_b2_all_nonimmune)
FeaturePlot_scCustom(b1_b2_all_nonimmune, features = "Pecam1", split.by = "Age")
FeaturePlot_scCustom(b1_b2_all_nonimmune, features = "Pecam1")

Idents(b1_b2_all_nonimmune) <- "Age"
p5_b1_b2_all_nonimmune <- subset(b1_b2_all_nonimmune, idents = "P5")

p5_b1_b2_all_nonimmune <- NormalizeData(object = p5_b1_b2_all_nonimmune)
p5_b1_b2_all_nonimmune <- FindVariableFeatures(object = p5_b1_b2_all_nonimmune)
p5_b1_b2_all_nonimmune <- ScaleData(object = p5_b1_b2_all_nonimmune)

###PCA
p5_b1_b2_all_nonimmune <- RunPCA(object = p5_b1_b2_all_nonimmune)
print(x = p5_b1_b2_all_nonimmune[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = p5_b1_b2_all_nonimmune, dims = 1:6)

DimPlot(object = p5_b1_b2_all_nonimmune)
p5_b1_b2_all_nonimmune <- ProjectDim(object = p5_b1_b2_all_nonimmune)

##JackStraw
p5_b1_b2_all_nonimmune <- JackStraw(object = p5_b1_b2_all_nonimmune, num.replicate = 100)
p5_b1_b2_all_nonimmune <- ScoreJackStraw(object = p5_b1_b2_all_nonimmune, dims = 1:20)
JackStrawPlot(object = p5_b1_b2_all_nonimmune, dims = 1:20)

ElbowPlot(object = p5_b1_b2_all_nonimmune, ndims = 20)

###
p5_b1_b2_all_nonimmune <- FindNeighbors(object = p5_b1_b2_all_nonimmune, dims = 1:9)
p5_b1_b2_all_nonimmune <- FindClusters(object = p5_b1_b2_all_nonimmune, resolution = 0.7)
##
p5_b1_b2_all_nonimmune <- RunUMAP(p5_b1_b2_all_nonimmune, dims = 1:9)

FeaturePlot(p5_b1_b2_all_nonimmune, features = "Col1a1", label = T, order = T)

Idents(p5_b1_b2_all_nonimmune) <- "seurat_clusters"
p5_fibros <- subset(p5_b1_b2_all_nonimmune, idents = c("0", "1", "2", "3", "5", "6", "8", "10"))

p5_fibros <- NormalizeData(object = p5_fibros)
p5_fibros <- FindVariableFeatures(object = p5_fibros)
p5_fibros <- ScaleData(object = p5_fibros)

###PCA
p5_fibros <- RunPCA(object = p5_fibros)
print(x = p5_fibros[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = p5_fibros, dims = 1:6)

DimPlot(object = p5_fibros)
p5_fibros <- ProjectDim(object = p5_fibros)

##JackStraw
p5_fibros <- JackStraw(object = p5_fibros, num.replicate = 100)
p5_fibros <- ScoreJackStraw(object = p5_fibros, dims = 1:20)
JackStrawPlot(object = p5_fibros, dims = 1:20)

ElbowPlot(object = p5_fibros, ndims = 20)

###
p5_fibros <- FindNeighbors(object = p5_fibros, dims = 1:9)
p5_fibros <- FindClusters(object = p5_fibros, resolution = 0.7)
##
p5_fibros <- RunUMAP(p5_fibros, dims = 1:9)

FeaturePlot(p5_fibros, features = "Col1a1", label = T, order = T)

p5_fibros$cell_type <- "Fibros"
qsave(p5_fibros, "/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/rds_files/p5_fibros_fibros_timecourse.rds")

###iTALK on adult ILCs, T cells, and Fibroblasts

p5_ILC_Tcells_fibros <- merge(x = p5_ILC_Tcells, y = p5_fibros)

p5_ILC_Tcells_fibros <- NormalizeData(object = p5_ILC_Tcells_fibros)
p5_ILC_Tcells_fibros <- FindVariableFeatures(object = p5_ILC_Tcells_fibros)
#length(x = VariableFeatures(object = all_immune))
p5_ILC_Tcells_fibros <- ScaleData(object = p5_ILC_Tcells_fibros)

###PCA
p5_ILC_Tcells_fibros <- RunPCA(object = p5_ILC_Tcells_fibros)
print(x = p5_ILC_Tcells_fibros[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = p5_ILC_Tcells_fibros, dims = 1:6)

DimPlot(object = p5_ILC_Tcells_fibros)
p5_ILC_Tcells_fibros <- ProjectDim(object = p5_ILC_Tcells_fibros)

##JackStraw
p5_ILC_Tcells_fibros <- JackStraw(object = p5_ILC_Tcells_fibros, num.replicate = 100)
p5_ILC_Tcells_fibros <- ScoreJackStraw(object = p5_ILC_Tcells_fibros, dims = 1:20)
JackStrawPlot(object = p5_ILC_Tcells_fibros, dims = 1:20)

ElbowPlot(object = p5_ILC_Tcells_fibros, ndims = 20)

###
p5_ILC_Tcells_fibros <- FindNeighbors(object = p5_ILC_Tcells_fibros, dims = 1:5)
p5_ILC_Tcells_fibros <- FindClusters(object = p5_ILC_Tcells_fibros, resolution = 0.7)
##
p5_ILC_Tcells_fibros <- RunUMAP(p5_ILC_Tcells_fibros, dims = 1:5)

qsave(p5_ILC_Tcells_fibros, "/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/rds_files/p5_ILC_Tcells_fibros_timecourse.rds")

Idents(p5_ILC_Tcells_fibros) <- "cell_type"

DimPlot_scCustom(p5_ILC_Tcells_fibros, label = T)

FeaturePlot_scCustom(p5_ILC_Tcells_fibros, features = "Gata3", label = T)

df <- t(as.matrix(GetAssayData(object = p5_ILC_Tcells_fibros, slot = "counts")))

colnames(df) <- toupper(colnames(df))


ct <- as.matrix(p5_ILC_Tcells_fibros$cell_type)

colnames(ct) <- "cell_type"

merged <- cbind(df, ct)
write.table(merged, '/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/P5_2geneXcell_Matrix_ILC_Tcells_fibros.txt', 
            sep = '\t', row.names = T, col.names = T, quote = F)

all_immune <- qread("/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/rds_files/timecourse_all_immune.rds")

Idents(b1_b2_all_nonimmune) <- "seurat_clusters"
DimPlot(b1_b2_all_nonimmune, split.by = "Age") + DarkTheme()

p5_ILC_TCell_fibros_iTALK <-read.table('/Users/tlanser/Desktop/Stevens Lab/Alec Single Cell Projects/immune_compartment_aging_project/data/2geneXcell_Matrix_ILC_Tcells_fibros_TC.txt', sep='\t', header=T, stringsAsFactors = F)
## highly expressed ligand-receptor pairs

# find top 50 percent highly expressed genes
p5_highly_exprs_genes<-rawParse(p5_ILC_TCell_fibros_iTALK,top_genes=50,stats='mean')
# find the ligand-receptor pairs from highly expressed genes
comm_list<-c('growth factor')
cell_col<-structure(c('#4a84ad','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b','#9ef49a'),names=unique(p5_ILC_TCell_fibros_iTALK$cell_type))
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

LRPlot(res_cat[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_mean_exprs[1:20],link.arr.width=res_cat$cell_to_mean_exprs[1:20])
title(comm_type)
res<-rbind(res,res_cat)
res<-res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:20,]
NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
LRPlot(res[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res$cell_from_mean_exprs[1:20],link.arr.width=res$cell_to_mean_exprs[1:20])

## significant ligand-receptor pairs between compare groups

# randomly assign the compare group to each sample
p5_ILC_TCell_fibros_iTALK<-p5_ILC_TCell_fibros_iTALK %>% mutate(compare_group=sample(2,nrow(p5_ILC_TCell_fibros_iTALK),replace=TRUE))
# find DEGenes of gdT cells and NK cells between these 2 groups
deg_t<-DEG(p5_ILC_TCell_fibros_iTALK %>% filter(cell_type=='gdTcells'),method='Wilcox',contrast=c(2,1))
deg_fibros<-DEG(p5_ILC_TCell_fibros_iTALK %>% filter(cell_type=='Fibros'),method='Wilcox',contrast=c(2,1))
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





