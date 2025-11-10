setwd("~/Documents/Immune")
library(Seurat)
library(SeuratDisk)
library(reticulate)
library(anndata)
library(rhdf5)
Convert("cell_com_umap.h5ad", dest = "h5seurat", overwrite = F)
cell_com<- LoadH5Seurat("cell_com_umap.h5seurat", assays="RNA", meta.data = FALSE, misc = FALSE)
obs <- h5read("cell_com_umap.h5ad", "/obs")
meta <- data.frame(lapply(names(obs), function(x) { 
  if (length(obs[[x]])==2) 
    obs[[x]][['categories']][ifelse(obs[[x]][['codes']] >= 0, obs[[x]][['codes']] + 1, NA)]
  else 
    as.numeric(obs[[x]])
}
), row.names=Cells(cell_com))#cambiare oggetto
colnames(meta) <- names(obs)

cell_com<- AddMetaData(cell_com,meta)
uns <- h5read("cell_com_umap.h5ad", "/uns") ####per avere dati gene markers

###celltype in new_celltype

#NICHENET'S
##analisi cell-cell interaction con oggetto Seurat
library(nichenetr)
library(tidyverse)
cell_com
#An object of class Seurat 
#27836 features across 15615 samples within 1 assay 
#Active assay: RNA (27836 features, 0 variable features)
#3 layers present: counts, data, scale.data
#16 dimensional reductions calculated: 

DefaultAssay(cell_com)<- "RNA"
DimPlot(cell_com, reduction = "umap", group.by = "new_celltype")
cell_com@meta.data$new_celltype %>% table()
#Anxa       Bcell         CNT       CNT_2        Cyst    INF_resp     KO_GSTO      Macs_1   Macs_Mrc1     Mono_DC 
#408         173         513         144        2321         139        2983        3699         669         513 
#Mono_DC2          Na Neutrophils          PT        T_NK       Tcell     Uretric 
#422          20        1433         607        1256          48         267 

cell_com<- alias_to_symbol_seurat(cell_com, "mouse")
Idents(cell_com)<- "new_celltype"
#Read in NicheNet's ligand-target prior model, ligand-receptor network and weighted integrated networks:

organism <- "mouse"

if(organism == "human"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
} else if(organism == "mouse"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
  
}

lr_network <- lr_network %>% distinct(from, to)
head(lr_network)

#ricordarsi di impostare come assay RNA
# seleziono i geni espressi almeno al 5% in un cluster. 
receiver = "Macs_Mrc1"
expressed_genes_receiver <- get_expressed_genes(receiver, cell_com, pct = 0.05)

#Get a list of all receptors available in the ligand-receptor network, 
#and define expressed receptors as genes that are in the ligand-receptor network and expressed in the receiver. 
#Then, define the potential ligands as all ligands whose cognate receptors are expressed.

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

sender_celltypes <- c("Cyst")

# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender= uns_matrix$names.10[uns_matrix$pvals_adj.5<= 0.05 & uns_matrix$logfoldchanges.5 >= 1] 

list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, cell_com, 0.1)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# Also check 
length(expressed_genes_sender)
## [1] 9948
length(potential_ligands)
## [1] 750
length(potential_ligands_focused)
## [1] 163
VlnPlot(cell_com, features= potential_ligands_focused[1:10])
#Define the gene set of interest
##solitamente DEG in questo caso uso marker genes
###
###
uns_matrix<- as.data.frame(uns$rank_genes_groups)

marker= uns_matrix$names.5[uns_matrix$pvals_adj.5<= 0.05 & uns_matrix$logfoldchanges.5 >= 1] 
length(marker)
geneset_oi = marker %>% .[. %in% rownames(ligand_target_matrix)]
length(geneset_oi)
#[1] 440
#3. Define the background genes

background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
length(background_expressed_genes)
#[1] 6987

#4. Perform NicheNet ligand activity analysis
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands_focused)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities
p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

head(active_ligand_target_links_df)

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

