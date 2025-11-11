.libPaths("/home/emilia/miniconda3/envs/scrna/lib/R/library")
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
setwd("/mnt/storage/scrna/data/")
# Load the datasets
Lib_sampleno85274_data <- Read10X(data.dir =  "/mnt/storage/scrna/data/Lib_sampleno85274/outs/filtered_feature_bc_matrix")
Lib_sampleno85275_data <- Read10X(data.dir =  "/mnt/storage/scrna/data/Lib_sampleno85275/outs/filtered_feature_bc_matrix")
Lib_sampleno85276_data <- Read10X(data.dir =  "/mnt/storage/scrna/data/Lib_sampleno85276/outs/filtered_feature_bc_matrix")
Lib_sampleno85277_data <- Read10X(data.dir =  "/mnt/storage/scrna/data/Lib_sampleno85277/outs/filtered_feature_bc_matrix")
Lib_sampleno85278_data <- Read10X(data.dir =  "/mnt/storage/scrna/data/Lib_sampleno85278/outs/filtered_feature_bc_matrix")
Lib_sampleno85279_data <- Read10X(data.dir =  "/mnt/storage/scrna/data/Lib_sampleno85279/outs/filtered_feature_bc_matrix")



# Initialize the Seurat object with the raw (non-normalized data).
Lib_sampleno85274 <- CreateSeuratObject(counts = Lib_sampleno85274_data, project = "Lib_sampleno85274", min.cells = 3, min.features = 200)
Lib_sampleno85275 <- CreateSeuratObject(counts = Lib_sampleno85275_data, project = "Lib_sampleno85275", min.cells = 3, min.features = 200)
Lib_sampleno85276 <- CreateSeuratObject(counts = Lib_sampleno85276_data, project = "Lib_sampleno85276", min.cells = 3, min.features = 200)
Lib_sampleno85277 <- CreateSeuratObject(counts = Lib_sampleno85277_data, project = "Lib_sampleno85277", min.cells = 3, min.features = 200)
Lib_sampleno85278 <- CreateSeuratObject(counts = Lib_sampleno85278_data, project = "Lib_sampleno85278", min.cells = 3, min.features = 200)
Lib_sampleno85279 <- CreateSeuratObject(counts = Lib_sampleno85279_data, project = "Lib_sampleno85279", min.cells = 3, min.features = 200)


Lib_sampleno8527<- merge(Lib_sampleno85274, y = c(Lib_sampleno85275,Lib_sampleno85276,Lib_sampleno85277,Lib_sampleno85278,Lib_sampleno85279), add.cell.ids = c("MZ74","MZ75","MZ76","MZ77","MZ78","MZ79"), project = "Lib_sampleno8527")
#Lib_sampleno8527@meta.data$dataset <- c(rep("Lib_sampleno85274", ncol(Lib_sampleno85274)), rep("Lib_sampleno85275", ncol(Lib_sampleno85275)),rep("Lib_sampleno85276", ncol(Lib_sampleno85276)),rep("Lib_sampleno85277", ncol(Lib_sampleno85277)),rep("Lib_sampleno85278", ncol(Lib_sampleno85278)),rep("Lib_sampleno85279", ncol(Lib_sampleno85279)))


Lib_sampleno8527[["percent.mt"]] <- PercentageFeatureSet(Lib_sampleno8527, pattern = "^mt-")

dir.create("./results_seurat_harmony/test")
setwd("./results_seurat_harmony/test")
pdf("vln.pdf")
VlnPlot(Lib_sampleno8527, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


pdf("QCs.pdf")
plot1 <- FeatureScatter(Lib_sampleno8527, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Lib_sampleno8527, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()


Lib_sampleno8527 <- subset(Lib_sampleno8527, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 1)



Lib_sampleno8527 <- Lib_sampleno8527 %>% Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = Lib_sampleno8527@var.genes, npcs = 20, verbose = FALSE)



pdf("pcs_before_harmony.pdf",width=811/72, height=478/72)    
p1 <- DimPlot(object = Lib_sampleno8527, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- DimPlot(object = Lib_sampleno8527, reduction = "pca", pt.size = .1, group.by = "orig.ident", split.by="orig.ident")
p3<- VlnPlot(object = Lib_sampleno8527, features = "PC_1", group.by = "orig.ident",  pt.size = .1)
p1
p2
p3
dev.off()

Lib_sampleno8527_js <- JackStraw(Lib_sampleno8527, num.replicate = 100)
Lib_sampleno8527_js <- ScoreJackStraw(Lib_sampleno8527_js, dims = 1:20)


pdf("JackStrawPlot.pdf")
JackStrawPlot(Lib_sampleno8527_js, dims = 1:20)
dev.off()



library(harmony)



Lib_sampleno8527_harmony <- Lib_sampleno8527 %>% RunHarmony( "orig.ident", plot_convergence = F)

pdf("pcs_after_harmony.pdf",width=811/72, height=478/72)    
p1 <- DimPlot(object = Lib_sampleno8527_harmony, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- DimPlot(object = Lib_sampleno8527_harmony, reduction = "harmony", pt.size = .1, group.by = "orig.ident", split.by="orig.ident")
p3 <- VlnPlot(object = Lib_sampleno8527_harmony, features = "harmony_1", group.by = "orig.ident",  pt.size = .1)
p1
p2
p3
dev.off()


########################## Clustree ##################################################
Lib_sampleno8527_clust <- FindNeighbors(Lib_sampleno8527_harmony, reduction = 'harmony', dims = 1:20) %>% FindClusters(resolution = seq(from=0, to=1, by=0.1), print.output = 0, save.SNN = F)

pdf('Resolution_tree_up_to_1_pca_1_20.pdf',  width=17, height=17)
clustree(Lib_sampleno8527_clust)
dev.off()


######################################################################



Lib_sampleno8527_harmony <- Lib_sampleno8527_harmony %>%  
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.1) %>% 
    identity()

freq_table<-table(Lib_sampleno8527_harmony@meta.data$RNA_snn_res.0.1, Lib_sampleno8527_harmony@meta.data$orig.ident)
write.table(freq_table, "freq_table.txt", sep="\t")


pdf("feature_plot.pdf",width=811/72, height=478/72)    
FeaturePlot(Lib_sampleno8527_harmony, features = c("Foxp3", "Itgae", "Cd4", "Cd8a"))
dev.off()


pdf("VlnPlot_genes.pdf")
VlnPlot(Lib_sampleno8527_harmony, features = c("Foxp3", "Itgae", "Cd4", "Cd8a"))
dev.off()

pdf("dotplot.pdf")
DotPlot(Lib_sampleno8527_harmony, features = c("Foxp3", "Itgae", "Cd4", "Cd8a")) + RotatedAxis()
dev.off()



pdf("plot_umap_after_harmony.pdf",width=811/72, height=478/72)    

DimPlot(Lib_sampleno8527_harmony, reduction = "umap", group.by = "orig.ident", pt.size = .1)
DimPlot(Lib_sampleno8527_harmony, reduction = "umap", pt.size = .1, group.by = "orig.ident", split.by="orig.ident")

dev.off()

Lib_sampleno8527<-RunUMAP(Lib_sampleno8527, dims = 1:20)    
pdf("plot_umap_before_harmony.pdf",width=811/72, height=478/72)    

DimPlot(Lib_sampleno8527, reduction = "umap", group.by = "orig.ident", pt.size = .1)
DimPlot(Lib_sampleno8527, reduction = "umap", group.by = "orig.ident", pt.size = .1,split.by="orig.ident")

dev.off()

save(list=ls(),file="ws_merge_harmony.rda")









