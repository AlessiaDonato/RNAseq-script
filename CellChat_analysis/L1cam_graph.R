Idents(seurat_object)<- "celltype"
plot_density(seurat_object, "L1cam", reduction = "umap")
ggsave("density_L1cam_umap.png", width = 10, height = 8)
guitar_obj<- subset(seurat_object, idents = c("transforming_1", "transforming_2", 
                    "transforming_3", "transforming_4", "transforming_5",
                    "transforming_6", "transforming_7", "KO_Gsto1_1", "KO_Gsto1_2", 
                    "KO_Gsto1_3", "KO_IC", "KO_transitional_regular"))
plot_density(guitar_obj, "L1cam", reduction = "umap")
FeaturePlot(seurat_object, "L1cam", label =T)
ggsave("L1cam_expr.png",  width = 10, height = 8)
DotPlot(seurat_object, "L1cam") 
ggsave("L1cam_dotplot.png", width = 7, height = 10)

gene_name <- "L1cam"
gene_expr <- FetchData(seurat_object, vars = gene_name)
#gene_expr_2 <- FetchData(seurat_object, vars = gene_name)

# Plot densità
ggplot(gene_expr, aes_string(x = gene_name)) +
  geom_density(fill = "skyblue", alpha = 0.6, color = "black") +
  theme_minimal() +
  labs(title = paste(gene_name, "distribution"),
       x = "log-norm",
       y = "Density")+
  theme(
    axis.title = element_text(size = 16),    
    axis.text = element_text(size = 14),      
    plot.title = element_text(size = 18, face = "bold")  
  )
ggsave("L1cam_distr.png", width = 10, height = 8)

ggplot(gene_expr_2, aes_string(x = gene_name)) +
  geom_histogram(bins = 50, fill = "orange", color = "black") +
  theme_minimal() +
  labs(title = paste("Istogramma espressione di", gene_name),
       x = "Livello di espressione",
       y = "Numero di cellule")


# Separiamo i valori > 0
nonzero_values <- gene_expr$L1cam[gene_expr$L1cam > 0]

plot(density(nonzero_values),
     main = "Densità dei valori non zero di L1cam",
     xlab = "Valori di espressione",
     ylab = "Densità",
     col = "blue",
     lwd = 2)

# Calcoliamo i quantili per i valori > 0
q <- quantile(nonzero_values, probs = c(0.25, 0.75))
q
#     25%      75% 
#4.074101 5.034313

png("L1cam_nonzero.png", width = 800, height = 600)
hist(nonzero_values,
     breaks = 30,
     col = "lightblue",
     main = "L1cam non-zero distribution",
     xlab = "Expr",
     ylab = "Frequence", 
     cex.main = 1.5,  
     cex.axis = 1.5,  
     cex.lab = 1.5 )
abline(v = q, col = "red", lty = 2)
text(x = q, y = par("usr")[4], labels = names(q), pos = 3, cex = 0.8, col = "red")
dev.off()

groups <- rep("zero", length(gene_expr$L1cam))
groups[gene_expr$L1cam > 0 & gene_expr$L1cam <= q[1]] <- "low"
groups[gene_expr$L1cam > q[1] & gene_expr$L1cam <= q[2]] <- "medium"
#groups[gene_expr$L1cam > q[2] & gene_expr$L1cam <= q[3]] <- "media-alta"
groups[gene_expr$L1cam > q[2]] <- "high"
# Controllo la distribuzione
groups <- factor(groups, levels = c("zero", "low", "medium", "high"))
table(groups)
#groups
#zero    low medium   high 
#51758   2786   5570   2786 

seurat_object$L1cam_group <- groups
DimPlot(seurat_object, reduction = "umap", group.by = "L1cam_group")
FeaturePlot(seurat_object, features = "L1cam") + 
  scale_color_gradientn(colors = c("grey", "blue", "red"))

library(ggplot2)

umap_df <- Embeddings(seurat_object, "umap") %>% as.data.frame()
umap_df$group <- seurat_object$L1cam_group

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = group)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c("zero" = "grey", "low" = "blue", 
                                "medium" = "orange",
                                "high" = "red")) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 12),    # dimensione testo legenda
    legend.title = element_text(size = 14)    # dimensione titolo legenda, se presente
  )
ggsave("l1cam_quantile.png", width = 10, height = 8)
DimPlot(seurat_object, reduction = "umap", group.by = "celltype",label = T) +
  NoLegend()
ggsave("umap_celltype.png", width = 10, height = 8)


mu <- mean(nonzero_values)
sigma <- sd(nonzero_values)
groups <- rep("zero", length(gene_expr$L1cam))

groups[gene_expr$L1cam > 0 & gene_expr$L1cam < (mu - 2*sigma)] <- "bassa"
groups[gene_expr$L1cam > (mu - 2*sigma) & gene_expr$L1cam < (mu + 2*sigma)] <- "media"

groups[gene_expr$L1cam > (mu + 2*sigma)] <- "alta"

# Controllo
table(groups)

hist(nonzero_values, breaks = 40, col = "lightblue", main = "Distribuzione non-zero",
     xlab = "Espressione")
abline(v = c(mu + 2*sigma,mu - 2*sigma) , col = "red", lwd = 2, lty = 2)
legend("topright", legend = "media + 2sd", col = "red", lty = 2, lwd = 2)

# 90% CI
groups <- rep("zero", length(gene_expr$L1cam))  

groups[gene_expr$L1cam > 0 & gene_expr$L1cam < (mu - 1.645 * sigma)] <- "low"
groups[gene_expr$L1cam > (mu - 1.645 * sigma) & gene_expr$L1cam < (mu + 1.645*sigma)] <- "medium"
groups[gene_expr$L1cam > (mu + 1.645 * sigma)] <- "high"

# Controllo
table(groups)
#groups
#high    low medium   zero 
# 633    437  10072  51758

hist(nonzero_values, breaks = 40, col = "lightblue", main = "Distribuzione non-zero",
     xlab = "Espressione")
abline(v = c(mu + 1.645*sigma,mu - 1.645*sigma) , col = "red", lwd = 2, lty = 2)
legend("topright", legend = "media + 2sd", col = "red", lty = 2, lwd = 2)


seurat_object$L1cam_group_CI90 <- groups

umap_df <- Embeddings(seurat_object, "umap") %>% as.data.frame()
umap_df$group <- seurat_object$L1cam_group_CI90

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = group)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c("zero" = "grey", "low" = "blue", 
                                "medium" = "orange",
                                "high" = "red")) +
  theme_minimal() 
ggsave("l1cam_quantile.png", width = 10, height = 8)

######Nt5e
plot_density(guitar_obj, "Nt5e", reduction = "umap")
FeaturePlot(guitar_obj, features = "Nt5e")
gene_Nt5e <- "Nt5e"
gene_expr <- FetchData(seurat_object, vars = gene_Nt5e)
#gene_expr_2 <- FetchData(seurat_object, vars = gene_name)

# Plot densità
ggplot(gene_expr, aes_string(x = gene_Nt5e)) +
  geom_density(fill = "skyblue", alpha = 0.6, color = "black") +
  theme_minimal() +
  labs(title = paste(gene_Nt5e, "distribution"),
       x = "log-norm",
       y = "Density")+
  theme(
    axis.title = element_text(size = 16),    
    axis.text = element_text(size = 14),      
    plot.title = element_text(size = 18, face = "bold")  
  )
ggsave("Nt5e_distr.png", width = 10, height = 8)

ggplot(gene_expr, aes_string(x = gene_Nt5e)) +
  geom_histogram(bins = 50, fill = "orange", color = "black") +
  theme_minimal() +
  labs(title = paste("Istogramma espressione di", gene_name),
       x = "Livello di espressione",
       y = "Numero di cellule")


# Separiamo i valori > 0
nonzero_values <- gene_expr$Nt5e[gene_expr$Nt5e > 0]

plot(density(nonzero_values),
     main = "Densità dei valori non zero di L1cam",
     xlab = "Valori di espressione",
     ylab = "Densità",
     col = "blue",
     lwd = 2)

# Calcoliamo i quantili per i valori > 0
q <- quantile(nonzero_values, probs = c(0.25, 0.75))
q
#     25%      75% 
#4.082977 5.135206

png("L1cam_nonzero.png", width = 800, height = 600)
hist(nonzero_values,
     breaks = 30,
     col = "lightblue",
     main = "L1cam non-zero distribution",
     xlab = "Expr",
     ylab = "Frequence", 
     cex.main = 1.5,  
     cex.axis = 1.5,  
     cex.lab = 1.5 )
abline(v = q, col = "red", lty = 2)
text(x = q, y = par("usr")[4], labels = names(q), pos = 3, cex = 0.8, col = "red")
dev.off()

groups <- rep("zero", length(gene_expr$Nt5e))
groups[gene_expr$Nt5e > 0 & gene_expr$Nt5e <= q[1]] <- "low"
groups[gene_expr$Nt5e > q[1] & gene_expr$Nt5e <= q[2]] <- "medium"
#groups[gene_expr$L1cam > q[2] & gene_expr$L1cam <= q[3]] <- "media-alta"
groups[gene_expr$Nt5e > q[2]] <- "high"
# Controllo la distribuzione
groups <- factor(groups, levels = c("zero", "low", "medium", "high"))
table(groups)
#groups
#zero    low medium   high 
#50238   3166   6330   3166  

seurat_object$Nt5e_group <- groups
DimPlot(seurat_object, reduction = "umap", group.by = "Nt5e_group")
FeaturePlot(seurat_object, features = "Nt5e") + 
  scale_color_gradientn(colors = c("grey", "blue", "red"))

library(ggplot2)

umap_df <- Embeddings(seurat_object, "umap") %>% as.data.frame()
umap_df$group <- seurat_object$Nt5e_group

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = group)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c("zero" = "grey", "low" = "blue", 
                                "medium" = "orange",
                                "high" = "red")) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 12),    # dimensione testo legenda
    legend.title = element_text(size = 14)    # dimensione titolo legenda, se presente
  )
ggsave("l1cam_quantile.png", width = 10, height = 8)
DimPlot(seurat_object, reduction = "umap", group.by = "celltype",label = T) +
  NoLegend()
ggsave("umap_celltype.png", width = 10, height = 8)

####
expr <- FetchData(guitar_obj, vars = c("L1cam", "Nt5e"))
head(expr)
dim(expr)

ggplot(expr, aes(x = L1cam, y = Nt5e)) +
  geom_point(alpha = 0.4, size = 0.6) +
  theme_minimal() +
  labs(title = "Co-expression of GeneA and GeneB",
       x = "GeneA expression",
       y = "GeneB expression")

library(ggplot2)
library(ggpointdensity)  # Per aggiungere densità come colore

# Estrai espressione
expr <- FetchData(guitar_obj, vars = c("L1cam", "Nt5e"))

# Scatter plot con densità
ggplot(expr, aes(x = L1cam, y = Nt5e)) +
  geom_pointdensity(adjust = 1, size = 0.5) +
  scale_color_viridis_c(option = "D") +
  theme_minimal() +
  labs(title = "Co-expression of L1cam and Nt5e",
       x = "L1cam expression (log-normalized)",
       y = "Nt5e expression (log-normalized)",
       color = "Cell density")

expr$group <- "none"
expr$group[expr$L1cam > 0 & expr$Nt5e == 0] <- "L1cam only"
expr$group[expr$L1cam == 0 & expr$Nt5e > 0] <- "Nt5e only"
expr$group[expr$L1cam > 0 & expr$Nt5e > 0] <- "both"
table(expr$group)
#both L1cam only       none  Nt5e only 
#1344       3015       7719       2163

# Aggiungi al Seurat object
guitar_obj$coexpr_group <- expr$group

DimPlot(seurat_object, group.by = "coexpr_group", reduction = "umap", 
        cols = c("none" = "lightgrey", "L1cam only" = "blue", 
                 "Nt5e only" = "green", "both" = "red"))

umap_df <- Embeddings(guitar_obj, "umap") %>% as.data.frame()
umap_df$group <- guitar_obj$coexpr_group

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = group)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c("none" = "grey", "L1cam only" = "blue", 
                                "Nt5e only" = "green",
                                "both" = "red")) +
  theme_minimal() 

# Estrai espressione
expr <- FetchData(seurat_object, vars = c("L1cam", "Nt5e"))

# Scatter plot con densità
ggplot(expr, aes(x = L1cam, y = Nt5e)) +
  geom_pointdensity(adjust = 1, size = 0.5) +
  scale_color_viridis_c(option = "D") +
  theme_minimal() +
  labs(title = "Co-expression of L1cam and Nt5e",
       x = "L1cam expression (log-normalized)",
       y = "Nt5e expression (log-normalized)",
       color = "Cell density")

expr$group <- "none"
expr$group[expr$L1cam > 0 & expr$Nt5e == 0] <- "L1cam only"
expr$group[expr$L1cam == 0 & expr$Nt5e > 0] <- "Nt5e only"
expr$group[expr$L1cam > 0 & expr$Nt5e > 0] <- "both"
table(expr$group)

# Aggiungi al Seurat object
seurat_object$coexpr_group <- expr$group

DimPlot(seurat_object, group.by = "coexpr_group", reduction = "umap", 
        cols = c("none" = "lightgrey", "L1cam only" = "blue", 
                 "Nt5e only" = "green", "both" = "red"))


