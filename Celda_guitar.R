library(zellkonverter)
library(SingleCellExperiment)
library(celda)
library(scater)
all_cell <- readH5AD("~/Documents/Cellchat/harmony_P33_P80_integrated_hvg_3000_no_outliers_n10.h5ad")
dim(all_cell)
#[1] 27836 62900
levels(all_cell$celltype)
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
#[1] 27836 23347
colData(simsce)<- colData(simsce)[, 1:30]
table(colData(simsce)$celltype) 
library(RColorBrewer)
n <- length(unique(simsce$celltype)) 
cols <- colorRampPalette(brewer.pal(12, "Set3"))(n)
plotReducedDim(simsce, dimred = "X_umap", colour_by = "celltype") +
  scale_color_manual(values = cols)

simsce <- selectFeatures(simsce)
altExp(simsce, "featureSubset")
#dim: 11614 23347 

#to select number of feature module (L)
#set interval 
moduleSplit <- recursiveSplitModule(simsce, altExpName = "featureSubset", 
                                    initialL = 2, maxL = 15)
#elbowplot (Visualize perplexity of every model in a celdaList, by unique K/L combinations)
plotGridSearchPerplexity(moduleSplit)
plotRPC(moduleSplit)
#to select cell cluster number (K)
moduleSplitSelect <- subsetCeldaList(moduleSplit, params = list(L = 7))
cellSplit <- recursiveSplitCell(moduleSplitSelect,
                                initialK = 3,
                                maxK = 20,
                                yInit = celdaModules(moduleSplitSelect))
plotGridSearchPerplexity(cellSplit)
plotRPC(cellSplit)
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
my_colors <- c("1" = "red", "2" = "cyan", "3" = "orange", "4" = "blue", 
               "5" = "yellow")
png("celda_umap.png", width = 1000, height = 1000)
plotReducedDim(sce, dimred = "X_umap", colour_by = "celda_cluster") +
  scale_color_manual(values = my_colors)
dev.off()

feat_sub <- altExp(sce, "featureSubset")  # subset di feature usate
modules<- celdaModules(sce)
names(modules) <- rownames(feat_sub)
head(names(modules))
head(modules)
tab <- table(sce$celltype, sce$celda_cluster)
tab
prop_tab <- prop.table(tab, margin = 1) * 100
print(round(prop_tab, 2))
# Supponiamo che i cluster da tenere siano in un vettore
clusters_to_keep <- c("KO_Gsto1_1", "KO_Gsto1_2", "KO_Gsto1_3", 
                      "KO_IC", "KO_transitional_regular", 
                      "transforming_1", "transforming_2", 
                      "transforming_3", "transforming_4", 
                      "transforming_5", "transforming_6", 
                      "transforming_7","transforming_mixed",
                      "KO_CNT", "KO_DCT_1","KO_transitional_transforming",
                      "KO_DCT_2", "KO_PC_1", "KO_PC_2")

# Filtra righe
prop_tab_filtered <- prop_tab[rownames(prop_tab) %in% clusters_to_keep, ]
print(prop_tab_filtered)

library(ggalluvial)

df <- as.data.frame(as.table(prop_tab_filtered))
colnames(df) <- c("Old", "New", "Freq")

colnames(df) <- c("OldCluster", "NewCluster", "Freq")

# Verifica se è alluviale (serve per geom_alluvium)
is_alluvia_form(df, axes = 1:2, silent = TRUE)  # dovrebbe dare TRUE

# Crea il plot
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
# Se la tabella è già normalizzata (proporzioni)
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
   # Crea file PNG ad alta risoluzione
  png(paste0("module_heatmap_", mod, ".png"), width = 2000, height = 2000, 
      res = 300)
  
  # Crea e stampa la heatmap
  plot(moduleHeatmap(sce, featureModule = mod))
  
  # Chiude il file
  dev.off()
}

genes_mod2 <- names(modules)[modules == 2] 

#########################################

#simsce=readRDS("guitar_completo.rds")
guitar<- readRDS("guitar_completo.rds")

simsce <- selectFeatures(guitar, minCount = 5, minCell =50)
altExp(simsce, "featureSubset")
#dim: 5134 23347 

moduleSplit <- recursiveSplitModule(simsce, altExpName = "featureSubset", 
                                    initialL = 2, maxL = 15)
#elbowplot (Visualize perplexity of every model in a celdaList, by unique K/L combinations)
plotGridSearchPerplexity(moduleSplit, altExpName = "featureSubset")
plotRPC(moduleSplit, altExpName = "featureSubset")
#per selezionare numero moduli di cellule K
moduleSplitSelect <- subsetCeldaList(moduleSplit, altExpName = "featureSubset",
                                     params = list(L = 7))
cellSplit <- recursiveSplitCell(moduleSplitSelect,altExpName = "featureSubset",
                                initialK = 3,
                                maxK = 15,
                                yInit = celdaModules(moduleSplitSelect))
plotGridSearchPerplexity(cellSplit, altExpName = "featureSubset")
plotRPC(cellSplit,altExpName = "featureSubset")
sce <- celda_CG(x = simsce,altExpName ="featureSubset",  K = 5, L = 7, 
                verbose = FALSE, nchains = 5)

table(celdaClusters(sce), simsce$celltype)
colData(sce)$celda_cluster <- celdaClusters(sce)
plot(celdaHeatmap(sce = sce, nfeatures = 10, showNamesFeature=T, ))
my_colors <- c("1" = "red", "2" = "cyan", "3" = "orange", "4" = "blue", 
               "5" = "yellow")
#png("celda_umap.png", width = 1000, height = 1000)
plotReducedDim(sce, dimred = "X_umap", colour_by = "celda_cluster") +
  scale_color_manual(values = my_colors)
dev.off()

feat_sub <- altExp(sce, "featureSubset")  # subset di feature usate
modules<- celdaModules(sce)
names(modules) <- rownames(feat_sub)
head(names(modules))
head(modules)
tab <- table(sce$celltype, sce$celda_cluster)
tab
prop_tab <- prop.table(tab, margin = 1) * 100
print(round(prop_tab, 2))
# Supponiamo che i cluster da tenere siano in un vettore
clusters_to_keep <- c("KO_Gsto1_1", "KO_Gsto1_2", "KO_Gsto1_3", 
                      "KO_IC", "KO_transitional_regular", 
                      "transforming_1", "transforming_2", 
                      "transforming_3", "transforming_4", 
                      "transforming_5", "transforming_6", 
                      "transforming_7","transforming_mixed",
                      "KO_CNT", "KO_DCT_1","KO_transitional_transforming",
                      "KO_DCT_2", "KO_PC_1", "KO_PC_2")

# Filtra righe
prop_tab_filtered <- prop_tab[rownames(prop_tab) %in% clusters_to_keep, ]
print(prop_tab_filtered)

library(ggalluvial)

df <- as.data.frame(as.table(prop_tab_filtered))
colnames(df) <- c("Old", "New", "Freq")

colnames(df) <- c("OldCluster", "NewCluster", "Freq")

# Verifica se è alluviale (serve per geom_alluvium)
is_alluvia_form(df, axes = 1:2, silent = TRUE)  # dovrebbe dare TRUE

# Crea il plot
ggplot(df,
       aes(axis1 = OldCluster, axis2 = NewCluster, y = Freq)) +
  geom_alluvium(aes(fill = OldCluster), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("OldCluster", "NewCluster"),expand = c(.05, .05))+
  theme_minimal() +
  labs(title = "Flusso da vecchi cluster a nuovi (celda)",
       y = "Numero cellule")
ggsave("celda_prob_alluvial.png", width =12, height = 10)
dev.off()
# Se la tabella è già normalizzata (proporzioni)
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
  # Crea file PNG ad alta risoluzione
  png(paste0("module_heatmap_", mod, ".png"), width = 2000, height = 2000, 
      res = 300)
  
  # Crea e stampa la heatmap
  print(moduleHeatmap(sce, featureModule = mod))
  
  dev.off()
}

table(celdaModules(sce))
#  1    2    3    4    5    6    7 
#257  156  686 2470 1349   21  195

summary(celdaClusters(sce))
#   1    2    3    4    5 
#3704 5003 7012 4102 3526 

for (mod in unique_modules) {
  genes_mod <- names(modules)[modules == mod] 
  write.table(genes_mod,paste0("module_", mod, "_50counts_5cells.txt"))
}

saveRDS(sce, "celda_5counts_50cell.rds")
