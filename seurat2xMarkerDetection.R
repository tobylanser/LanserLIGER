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
library(qs)
library(scran)


setwd("C:/Users/Stevens Lab/Desktop/LanserLIGER")

Sys.setenv('R_MAX_VSIZE'=memory.size(max=TRUE))

fibro_meninges <- readRDS(("fibro_meninges.rds"))
mouseSSfibro_seurat <- readRDS("Mouse_SS_Fibro_scaled.rds")

meninges = qread("P14_CB_all_CD45neg_arranged.qs")
mouseSSfibro = qread("Mouse_SS_Fibro_arranged.qs")

meninges=meninges[!is.na(rowData(meninges)$ensembl_gene_id),]
x=rowSums(counts(meninges))
meninges=meninges[order(x,decreasing = T),]
meninges=meninges[!duplicated(rowData(meninges)$ensembl_gene_id),]
row.names(meninges)=rowData(meninges)$ensembl_gene_id

nbt=Seurat::CreateSeuratObject(counts=counts(meninges),
                               project = "SeuratProject",
                               assay = "RNA",
                               min.cells = 0,
                               min.features = 0,
                               names.field = 1,
                               names.delim = "-",
                               meta.data = as.data.frame(colData(meninges)))

meninges_seurat <- nbt
meninges_seurat@meta.data$ClustName= "Dura"
table(meninges_seurat$ClustName,meninges_seurat$Tissue)

table(meninges_seurat$ClustName,meninges_seurat$Tissue)
table(meninges_seurat$Dura,meninges_seurat$Tissue)

mouseSSfibro=mouseSSfibro[!is.na(rowData(mouseSSfibro)$ensembl_gene_id),]
x=rowSums(counts(mouseSSfibro))
mouseSSfibro=mouseSSfibro[order(x,decreasing = T),]
mouseSSfibro=mouseSSfibro[!duplicated(rowData(mouseSSfibro)$ensembl_gene_id),]
row.names(mouseSSfibro)=rowData(mouseSSfibro)$ensembl_gene_id

nbt=Seurat::CreateSeuratObject(counts=counts(mouseSSfibro),
                               project = "SeuratProject",
                               assay = "RNA",
                               min.cells = 0,
                               min.features = 0,
                               names.field = 1,
                               names.delim = "-",
                               meta.data = as.data.frame(colData(mouseSSfibro)))

mouseSSfibro_seurat <- nbt
rm(nbt)


mouseSSfibro_seurat <- NormalizeData(mouseSSfibro_seurat)
mouseSSfibro_seurat <- ScaleData(mouseSSfibro_seurat)
saveRDS(mouseSSfibro_seurat, "Mouse_SS_Fibro_scaled.rds")

FeaturePlot(fibro_meninges, features = c("ENSMUSG00000032349"), raster = F)

# dim(meninges)
# class(meninges)
# head(rowData(meninges))
# head(colData(meninges))
# meninges_seurat= as.Seurat(meninges)
# meninges_seurat=.extraExport2SeuratFn(inputData=meninges,project = "scRNA")
# meninges_seurat= as.Seurat(meninges, counts = "counts")
# 
# nbt@assays$RNA@meta.features=as.data.frame(rowData(inputData = meninges))


mouseSSfibro_seurat$orig.ident <- "Turley"
meninges_seurat$orig.ident <- "Meninges"
# rm(nbt)

Idents(object = meninges_seurat) <- meninges_seurat@meta.data$"seurat_clusters"
fibro_meninges <- subset(meninges_seurat, idents = c("2", "9", "13", "14", "23", "24"))
fibro.2x <- merge(mouseSSfibro_seurat, y = fibro_meninges, project = "MergedFibro2x")



fibro.2x <- readRDS("fibro2x.rds")
Mouse_SS_Fibro_scaled.rds
fibro_meninges <- NormalizeData(fibro_meninges)
fibro_meninges <- ScaleData(fibro_meninges)

saveRDS(fibro_meninges, "fibro_meninges.rds")


# Set the identity of your cells to the desired column
Idents(object = fibro.2x) <- fibro.2x@meta.data$"Tissue"

fibro.2x_tissue_markers <- FindMarkers(fibro.2x, ident.1 = "Dura", group.by = "Tissue", only.pos = T)
fibro.2x_tissue_markers_200 <- fibro.2x_tissue_markers %>%
  top_n(200, avg_log2FC)
write.csv(fibro.2x_tissue_markers_200, "fibro.2x_tissue_markers_200.csv")


# Set the identity of your cells to the desired column
Idents(object = fibro.2x) <- fibro.2x@meta.data$"ClustName"

fibro.2x_clustName_markers <- FindMarkers(fibro.2x, ident.1 = "Dura", group.by = "ClustName", only.pos = T)
fibro.2x_clustName_markers_200 <- fibro.2x_clustName_markers %>%
  top_n(200, avg_log2FC)
write.csv(fibro.2x_clustName_markers_200, "fibro.2x_clustName_markers_200.csv")

table(meninges$ClustName,meninges$Tissue)

table(fibro_meninges$ClustName,fibro_meninges$Tissue)

fibro.2xMarkers <- fibro.2x
# set cutoff to zero expression
foxd1.pos <- subset(fibro.2xMarkers, subset = Foxd1 > 0)
foxd1.neg <- subset(fibro.2xMarkers, subset = Foxd1 = 0)
WhichCells(fibro.2xMarkers@assays[["RNA"]]@data["Foxd1",])>1
subset(fibro.2xMarkers, subset = Foxd1 > 0.32)
#create metadata foxd1 column
fodx1.pos.cells <- "Foxd1+"
fibro.2x@meta.data <- cbind(fibro.2x@meta.data, fodx1.pos.cells)
colnames(fibro.2x@meta.data)[which(names(fibro.2x@meta.data) == "fodx1.pos.cells")] <- "foxd1"

table(fibro.2x@meta.data$foxd1)
fibro.2xMarkers <- merge(foxd1.pos, y = foxd1.neg, project = "Foxd1Markers")
fibro.2x@meta.data$foxd1 <- fibro.2xMarkers@meta.data$foxd1

fodx1.neg.cells <- "Foxd1-"
foxd1.neg@meta.data <- cbind(foxd1.neg@meta.data, fodx1.neg.cells)
colnames(foxd1.neg@meta.data)[which(names(foxd1.neg@meta.data) == "fodx1.neg.cells")] <- "foxd1"

gene_ids_1 = c("ENSMUSG00000032349", "ENSMUSG00000050010", "ENSMUSG00000030351", "ENSMUSG00000037625",
             "ENSMUSG00000025776", "ENSMUSG00000027339", "ENSMUSG00000037049", "ENSMUSG00000039323", "ENSMUSG00000078302",
             "ENSMUSG00000033059")
gene_ids_2 =c("ENSMUSG00000033059", "ENSMUSG00000033208", "ENSMUSG00000028444",
             "ENSMUSG00000010307", "ENSMUSG00000032826", "ENSMUSG00000027340", "ENSMUSG00000032348", "ENSMUSG00000032717")
gene_ids_3 = c("ENSMUSG00000037656", "ENSMUSG00000026576", "ENSMUSG00000062691", "ENSMUSG00000029998", "ENSMUSG00000046714",
             "ENSMUSG00000022199", "ENSMUSG00000024597", "ENSMUSG00000087365", "ENSMUSG00000022091", "ENSMUSG00000026433")
gene_ids_4= c("ENSMUSG00000040794", "ENSMUSG00000001665", "ENSMUSG00000024665", "ENSMUSG00000020230", "ENSMUSG00000004415",
             "ENSMUSG00000052854", "ENSMUSG00000027217", "ENSMUSG00000037762", "ENSMUSG00000106874", "ENSMUSG00000041329")
gene_ids_5 = c("ENSMUSG00000041075", "ENSMUSG00000050471", "ENSMUSG00000040164", "ENSMUSG00000085399", "ENSMUSG00000025041",
             "ENSMUSG00000089661", "ENSMUSG00000035314", "ENSMUSG00000021719", "ENSMUSG00000056313", "ENSMUSG00000029821")
gene_ids_6 = c("ENSMUSG00000059146", "ENSMUSG00000026241", "ENSMUSG00000019851", "ENSMUSG00000030350", "ENSMUSG00000046593",
             "ENSMUSG00000026904", "ENSMUSG00000018900", "ENSMUSG00000035919", "ENSMUSG00000010122", "ENSMUSG00000071489",
             "ENSMUSG00000089669", "ENSMUSG00000027463", "ENSMUSG00000040055", "ENSMUSG00000060424", "ENSMUSG00000052188")

gene_names = c("Elovl5", "Shisa3", "Tspan11", "Cldn11", "Crispld1", "Rassf2", "Smpd1", "Igfbp2", "Foxd1", "Pygb")

fibro.2x@meta.data %>%
  group_by(Tissue,Phase) %>%
  count() %>%
  group_by(Tissue) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=Tissue,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cells per tissue expressing intersect genes")


avgExpress <- AverageExpression(fibro.2x, group.by = "Tissue", features = features)
AddModuleScore(fibro.2x, features = features)
DoHeatmap(object = mouseSSfibro_seurat, group.by = "Tissue", features = gene_ids_1)
Foxd1DotPlot <- DotPlot(object = fibro.2x, group.by = "Tissue", features = "ENSMUSG00000078302")
Foxd1DotPlot +  labs(title = "Foxd1 Expression Across Tissue")
VlnPlot(fibro.2x, group.by = "Tissue", features = "ENSMUSG00000078302")
gg_Fig$data

fibro_meninges$seurat_clusters=droplevels(fibro_meninges$seurat_clusters)

write.csv(gg_Fig$data, "pctExp_top10Intersect_seurat2X.csv")


PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}


calc_helper(fibro_meninges, genes = "ENSMUSG00000032349")
TissueExpression <- PrctCellExpringGene(fibro.2x, genes = "ENSMUSG00000056313", group.by = "Tissue")
barplot((100*TissueExpression$Cell_proportion), xlab = TissueExpression$Feature, ylab = "% Gene Expression", main ="Gene % Expression by Tissue", col = "dodgerblue3")

TissueExpression
table(TissueExpression$Cell_proportion)

foxd1Expression <- matrix(foxd1Expression)

fibro_meninges <- FindVariableFeatures(fibro_meninges)
fibro_meninges <- RunPCA(fibro_meninges)
fibro_meninges <- FindNeighbors(fibro_meninges, dims = 1:25)
fibro_meninges <- FindClusters(fibro_meninges, resolution = 1)
fibro_meninges <- RunUMAP(fibro_meninges, dims = 1:25)

install.packages(pheatmap)

# load package
library(pheatmap)

FeaturePlot(fibro_meninges, features = c("ENSMUSG00000032349"), raster = F, order = T, cols = c("black","white")) + DarkTheme()
pheatmap(fibro_meninges)

install.packages('animation')
library(animation)
FeaturePlot(fibro_meninges, features = "ENSMUSG00000050010", label = T,  split.by = 'Tissue')
fibro_meninges <- AddModuleScore(fibro_meninges, features = gene_ids, name = "gene_id_score")
FeaturePlot(fibro_meninges, features = gene_ids_1)
all_genes
dev.off()
DimPlot(fibro_meninges, group.by = "seurat_clusters")
DimPlot(fibro_meninges, label = T)
VlnPlot(mouseSSfibro_seurat, group.by ="Tissue", features = gene_ids_1)
DotPlot(seurat3x, group.by = "Tissue", features = "Foxd1", cols = c("black","white")) + RotatedAxis() + DarkTheme()

options(expressions=100000)

fibro_meninges_Markers <- FindAllMarkers(fibro_meninges, only.pos = T)
fibro_meninges_Markers_top_200 <- fibro_meninges_Markers %>% group_by(cluster) %>%
  top_n(200, avg_log2FC)
write.csv(fibro_meninges_Markers_top_200, "fibro_meninges_Markers_top_200.csv")


# Set the identity of your cells to the desired column
Idents(object = fibro_meninges) <- fibro_meninges@meta.data$"Tissue"

fibro_meninges_Markers_Tissue <- FindMarkers(fibro_meninges, ident.1 = "Dura", only.pos = T)
fibro_meninges_Markers_Tissue_top200 <- fibro_meninges_Markers_Tissue %>%
  top_n(200, avg_log2FC)
write.csv(fibro_meninges_Markers_Tissue_top200, "fibro_meninges_Markers_Tissue_top200.csv")



Tissue <- c("Bone", "Visceral Adipose", "Omentum", "Intestine", "Heart", "Subcutaneous Adipose",
            "Liver", "Lymph Node", "Lung", "Mesentery", "Muscle", "Pancreas", "Skin", "Spleen",
            "Tendon", "Artery", "Dura", "CalvBM", "TibBM")
Expression_Prct_Tissue <- c(Bone = 4.54046137, Visceral_Adipose = 0.12205045, Omntuum =0, Intestine = 0.02165674, Heart = 1.11015220, 
                      Subcutaneous_Adipose =0.07187781, Liver = 0, Lymph_Node = 0, Lung = 0.71796604,
                      Mesentery = 0, Muscle = 0.15772871, Pancreas = 0, Skin = 0.69134701, Spleen =0.01583531, Tendon =14.96598639,
                      Artery =2.87081340, Dura = 75.56753689,
                     CalvBM = 41.90476190, TibBM = 07.40740741)

barplot(Expression_Prct_Tissue ,xlab = "Tissue", ylab = "Foxd1 % Expression", main ="Foxd1 % Expression by Tissue", col = "dodgerblue3")
