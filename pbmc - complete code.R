##ScRNA seq pbmc cells

library(dplyr)
library(Seurat)
library(cowplot)
library(Matrix)
library(ggplot2)
library(sctransform)
getwd()

setwd("C:/Users/rraag/Downloads/Parent_NGSC3_DI_PBMC_filtered_feature_bc_matrix")
library1.data<-Read10X("filtered_feature_bc_matrix")
library1 <- CreateSeuratObject(counts = library1.data, min.cells = 10, min.features = 800, project = "library1")
library1$library <- "1"

library1
#######Quality conttol step#####################
Big_library <- PercentageFeatureSet(library1, pattern = "^mt-", col.name = "percent.mt")
head(Big_library@meta.data, 5)
tail(Big_library@meta.data, 5)
dim(Big_library)
#############
#######################

length(Big_library@meta.data$percent.mt)
VlnPlot(Big_library, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mean("nFeature_RNA") 
median(Big_library@meta.data$nFeature_RNA) 
plot1 <- FeatureScatter(Big_library, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Big_library, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

cowplot::plot_grid(plot1,plot2, labels = c("A","B")) 
ggsave("qualitycontrolcorrplots.tiff", plot = last_plot(), dpi = 800, width = 50, height = 15, units = "cm", device = "tiff") 
########################
summary(Big_library@meta.data$percent.mt)
############
#Subset retina object according to quality control metrics to remove unwanted cells from the dataset based on the number of features and percentage of mitochondrial counts 
mouse_dbt <- subset(Big_library, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5) 
dim(Big_library)
median(Big_library@meta.data$nFeature_RNA) 

######################Normalizing the data#######################3

Big_library <- NormalizeData(Big_library, normalization.method = "LogNormalize", scale.factor = 10000)
###Running line below will just overwrite normalization with default parameters
Big_library <- NormalizeData(Big_library)
#Access normalised values  
Big_library[["RNA"]]@data  

# run sctransform
Big_library <- FindVariableFeatures(Big_library, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Big_library), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Big_library)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

##############Scaling the data################################
all.genes <- rownames(Big_library)
all.genes
write.table(all.genes, "pbmc_genes.txt", sep="\t")
Big_library <- ScaleData(Big_library, vars.to.regress = "percent.mt")
#Access scaled data 
Big_library[["RNA"]]@scale.data[1:5,1:5] 
####################3linear dimensionality test#####################333
Big_library <- RunPCA(Big_library, features = VariableFeatures(object = Big_library))
print(Big_library[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Big_library, dims = 1:4, reduction = "pca")
DimPlot(Big_library, reduction = "pca")
DimHeatmap(Big_library, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Big_library, dims = 1:15, cells = 500, balanced = TRUE)
Big_library <- RunPCA(Big_library, npcs = 30, verbose = FALSE) #####samples
Big_library <- RunUMAP(Big_library, reduction = "pca", dims = 1:30)
p1 <- DimPlot(Big_library, reduction = "umap")
p2 <- DimPlot(Big_library, reduction = "umap", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2 ############

##########################dimensionality test############################

Big_library <- JackStraw(Big_library, num.replicate = 100, dims = 50)
Big_library <- ScoreJackStraw(Big_library, dims = 1:20)

JackStrawPlot(Big_library, dims = 1:15)
ElbowPlot(Big_library)

#############clustering cells##########################

Big_library <- FindNeighbors(Big_library, dims = 1:30)
Big_library <- FindClusters(Big_library, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(mouse_dbt), 5)
tail(Idents(Big_library), 10)
DimPlot(Big_library, label = TRUE) + NoLegend()

#######################UMAP/tSNE test###################33  
######UMAP/TSNE#######
Big_library <- RunUMAP(Big_library, dims = 1:30)
Big_library <- RunTSNE(Big_library, check_duplicates = FALSE, dims = 1:30)
DimPlot(Big_library, reduction = "tsne", label = TRUE)
DimPlot(Big_library, reduction = "umap", label = TRUE)
saveRDS(Big_library, file = "###New result for all pbmc tsne-umap result.rds")

#################clusters biomarkers#####################

###Need to set idents to cluster before finding markers
###Idents(object = Big_library) <- Big_library@meta.data$RNA_snn_res.0.5

cluster1.markers <- FindMarkers(Big_library, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 10)
cluster5.markers <- FindMarkers(Big_library, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 10)
cluster2.markers <- FindMarkers(Big_library, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 20)


Big_library_markers <- FindAllMarkers(Big_library, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Big_library_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
head(Big_library_markers,20)
write.table(Big_library_markers, "Big_library__markers pbmc output.txt", sep = "\t")
View(Big_library_markers)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
DimPlot(Big_library, reduction = "umap", label = TRUE)

FeaturePlot(Big_library, features = c("RBP7"), label = T)
FeaturePlot(Big_library, features = c("RPL5"), label = T)
FeaturePlot(Big_library, features = c("S100A4"), label = T)
FeaturePlot(Big_library, features = c("IL7R"), label = T)
FeaturePlot(Big_library, features = c("FCGR1A"), label = T)
FeaturePlot(Big_library, features = c("KLRD1"), label = T)
FeaturePlot(Big_library, features = c("FCGR3A"), label = T)
FeaturePlot(Big_library, features = c("MS4A1"), label = T)
FeaturePlot(Big_library, features = c("NKG7"), label = T)
FeaturePlot(Big_library, features = c("GNLY"), label = T)
FeaturePlot(Big_library, features = c("CD8B"), label = T)
FeaturePlot(Big_library, features = c("FCER1A"), label = T)
FeaturePlot(Big_library, features = c("RPS5"), label = T)

FeaturePlot(Big_library, features = c("S100A12"), label = T)
FeaturePlot(Big_library, features = c("TIGIT"), label = T)
FeaturePlot(Big_library, features = c("CST3"), label = T)
FeaturePlot(Big_library, features = c("CCR7"), label = T)
FeaturePlot(Big_library, features = c("IL7R"), label = T)
FeaturePlot(Big_library, features = c("LYZ"), label = T)
Big_library <- RenameIdents(Big_library,  '0'= "CD14+ Monocytes", '1' = " DC", '2' = "CD8+ T cell", '3' = "Naive CD4+ T", '4' = "Monocytic derived DC", 
                            '5' = "NK cell", '6' = "B cells", '7' = "FCGR3A+ Mono", '8' = "B cells", '9' = "Effector CD8+ memory T cell",'10' = "CD14+CD16+ Monocytes", '11' = "Treg cells", 
                            '12' = "Mono derived DC", '13' = "Naive CD4+ T", '14' = "Memory CD4+", '15' = "CD14+ Monocytes", '16' = "Naive CD4+ T")

DimPlot(Big_library, reduction = "umap", label = TRUE)

DimPlot(Big_library, reduction = "tsne", label = TRUE)

saveRDS(Big_library,file = "results_seuratpbmc_labelled.rds")   

levels(Big_library)

###number of cells
table(Idents(Big_library))
###proportions
prop.table(table(Idents(Big_library)))
