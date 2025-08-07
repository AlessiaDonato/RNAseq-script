# Load required libraries
library(anndata)
library(zellkonverter)
library(Seurat)
# Read the AnnData (.h5ad) file
adata<- read_h5ad("pathtofile.h5ad")
# Extract expression counts matrix and transpose it (cells as columns, genes as rows)
expr_data <- t(as.matrix(adata$layers["counts"]))
# Extract meta data
meta_data <- as.data.frame(adata$obs)[, 1:41] #select metadata to use

# Create Seurat object
seurat_object <- CreateSeuratObject(counts = expr_data,
                                  meta.data = meta_data)
#Seurat object summary
seurat_object 
# Remove temporary variables to free memory
rm(meta_data) 
rm(expr_data) 

# Extract scaled data from AnnData, transpose and assign to Seurat object
scaled_data <- as.matrix(adata$layers["scaled"])
scaled_data <- t(scaled_data)
seurat_object[["RNA"]]$scale.data<-scaled_data
rm(scaled_data)

# Extract log-normalized data from AnnData,transpose and assign to Seurat object
log_data <- as.matrix(adata$layers["data"])
log_data <- t(log_data)
seurat_object[["RNA"]]$data<-log_data
rm(log_data)

# Extract UMAP embeddings from AnnData
umap_data <- adata$obsm$X_umap
rownames(umap_data) <- rownames(adata$obs)
#Rename UMAP embedding columns; assumes 2D UMAP
colnames(umap_data) <- c("UMAP_1", "UMAP_2")  
#Create a DimReduc object for UMAP embeddings in Seurat
seurat_object[["umap"]]<- CreateDimReducObject(embeddings = umap_data, 
                                               key = "UMAP_", assay = "RNA")
# Extract PCA embeddings from AnnData
pca_data <- adata$obsm$X_pca
rownames(pca_data) <- rownames(adata$obs)

# Set column names for PCA components
# This works if there are exactly 50 components; otherwise, set dynamically 
# based on the number of columns
colnames(pca_data) <- c(1:50) 
seurat_object[["pca"]]<- CreateDimReducObject(embeddings = pca_data, 
                                               key = "PCA_", assay = "RNA")
## Clean up variables to free memory
rm(pca_data)
rm(umap_data)
