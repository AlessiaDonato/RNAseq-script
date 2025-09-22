library(zellkonverter)
library(SingleCellExperiment)
library(celda)
library(scater)

all_cell <- readH5AD("path_to_file")
dim(all_cell)

levels(all_cell$celltype)
#select clusters to use
simsce <-all_cell[,all_cell$celltype %in% c("KO_Gsto1_1", "KO_Gsto1_2", 
                                            "KO_Gsto1_3","KO_PC_1", "KO_PC_2",
                                            "KO_IC", "KO_transitional_regular", 
                                            "transforming_1", "transforming_2", 
                                            "transforming_3", "transforming_4", 
                                            "transforming_5", "transforming_6", 
                                            "transforming_7","transforming_mixed",
                                            "KO_CNT", "KO_DCT_1","KO_DCT_2",
                                            "KO_transitional_transforming")]

rm(all_cell)
dim(counts(simsce))
#select metadata to use
colData(simsce)<- colData(simsce)[, 1:30]
table(colData(simsce)$celltype) 

#to use to plot umap with high number of cluster
library(RColorBrewer)
#set the number of cluster and the palette to use
n <- length(unique(simsce$celltype)) 
cols <- colorRampPalette(brewer.pal(12, "Set3"))(n)
plotReducedDim(simsce, dimred = "X_umap", colour_by = "celltype") +
  scale_color_manual(values = cols)

#feature selection is performed to reduce the size of features used for clustering
#by default 3 counts in at least 3 cells are included
simsce <- selectFeatures(simsce)
altExp(simsce, "featureSubset")
#select 5 counts in at least 50 cells 
#simsce <- selectFeatures(guitar, minCount = 5, minCell =50)

# Explore optimal number of feature modules (L) by specifying a range
moduleSplit <- recursiveSplitModule(simsce, altExpName = "featureSubset", 
                                    initialL = 2, maxL = 15)

# Visualize perplexity across models for each K/L combination (elbow plot)
plotGridSearchPerplexity(moduleSplit)

# Plot rate of perplexity change (RPC) to aid in module selection
plotRPC(moduleSplit)

# Explore optimal number of cell clusters (K) by specifying a range
moduleSplitSelect <- subsetCeldaList(moduleSplit, params = list(L = 7))
cellSplit <- recursiveSplitCell(moduleSplitSelect,
                                initialK = 3,
                                maxK = 20,
                                yInit = celdaModules(moduleSplitSelect))
# Visualize perplexity across models for each K/L combination (elbow plot)
plotGridSearchPerplexity(cellSplit)

# Plot rate of perplexity change (RPC) to aid in module selection
plotRPC(cellSplit)

#celda_CG perform simultaneously cluster cells and features
#celda_C only cluster cells 
#celda_G only cluster features (gene modules)
sce <- celda_CG(x = simsce, K = 5, L = 7, verbose = FALSE, nchains = 5)

table(celdaClusters(sce), simsce$celltype)

sce$celda_cluster <- celdaClusters(sce)
reducedDimNames(sce)

library(scater)
plotReducedDim(sce, dimred = "X_umap", colour_by = "celda_cluster")
png("celda_heatmap.png", width = 1000, height = 1000)

plot(celdaHeatmap(sce = sce, nfeatures = 10, showNamesFeature=T))
dev.off()
colData(sce)$celda_cluster <- celdaClusters(sce)

# Generate a UMAP with the same cell cluster colors as the heatmap
my_colors <- c("1" = "red", "2" = "cyan", "3" = "orange", "4" = "blue", 
               "5" = "yellow")
png("celda_umap.png", width = 1000, height = 1000)
plotReducedDim(sce, dimred = "X_umap", colour_by = "celda_cluster") +
  scale_color_manual(values = my_colors)
dev.off()

#features subset 
feat_sub <- altExp(sce, "featureSubset")  
modules<- celdaModules(sce)
names(modules) <- rownames(feat_sub)
head(names(modules))
head(modules)
#tab with celda cluster and old celltype(cluster)
tab <- table(sce$celltype, sce$celda_cluster)
tab
prop_tab <- prop.table(tab, margin = 1) * 100
print(round(prop_tab, 2))
#if use a subset of cluster, specify the cluster for a filtered tab 
clusters_to_keep <- c("KO_Gsto1_1", "KO_Gsto1_2", "KO_Gsto1_3", 
                      "KO_IC", "KO_transitional_regular", 
                      "transforming_1", "transforming_2", 
                      "transforming_3", "transforming_4", 
                      "transforming_5", "transforming_6", 
                      "transforming_7","transforming_mixed",
                      "KO_CNT", "KO_DCT_1","KO_transitional_transforming",
                      "KO_DCT_2", "KO_PC_1", "KO_PC_2")

prop_tab_filtered <- prop_tab[rownames(prop_tab) %in% clusters_to_keep, ]
print(prop_tab_filtered)

#alluvial plot for old cluster and celda cluster 
library(ggalluvial)
df <- as.data.frame(as.table(prop_tab_filtered))
colnames(df) <- c("OldCluster", "NewCluster", "Freq")
# Check if the data frame is valid for alluvial plotting
is_alluvia_form(df, axes = 1:2, silent = TRUE)  

#
ggplot(df,
       aes(axis1 = OldCluster, axis2 = NewCluster, y = Freq)) +
  geom_alluvium(aes(fill = OldCluster), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("OldCluster","NewCluster"),expand = c(.05, .05)) +
  theme_minimal() +
  labs(title = "Flusso da vecchi cluster a nuovi (celda)",
       y = "Numero cellule")
ggsave("celda_prob_alluvial.png", width =12, height = 10)
dev.off()

# heatmap
heatmap_matrix <- round(prop.table(prop_tab_filtered, margin = 1) * 100, 1)

library(pheatmap)
png("heatmap_prob_celda.png", width = 1500, height = 2000, res = 300)
pheatmap(heatmap_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         main = "Distribuzione vecchi cluster nei nuovi (celda)")
dev.off()

moduleHeatmap(sce, featureModule = 7)
unique_modules <- sort(as.integer(levels(modules)))
for (mod in unique_modules) {
  png(paste0("module_heatmap_", mod, ".png"), width = 2000, height = 2000, 
      res = 300)
  plot(moduleHeatmap(sce, featureModule = mod))
  dev.off()
}

#genes in module
genes_mod2 <- names(modules)[modules == 2] 

for (mod in unique_modules) {
  genes_mod <- names(modules)[modules == mod] 
  write.table(genes_mod,paste0("module_", mod, "_50counts_5cells.txt"))
}

saveRDS(sce, "celda_5counts_50cell.rds")
