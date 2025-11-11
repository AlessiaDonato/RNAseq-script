library(CellChat)
library(patchwork)
library(ggplot2)

#####change pattern number
gene= intersect(rownames(cellchat@data), ligand_df$Pattern_3)
# filter gene
genes_present <- intersect(gene, rownames(cellchat@data))
genes_per_page <- 4
gene_chunks <- split(genes_present, ceiling(seq_along(genes_present) / genes_per_page))

#plot
for (i in seq_along(gene_chunks)) {
  chunk <- gene_chunks[[i]]
  plots <- lapply(chunk, function(g) {
    plotGeneExpression(cellchat, features = g, enriched.only = TRUE, type = "violin") +
      ggtitle(g)
  })
  
  if (length(plots) < 4) {
    plots <- c(plots, replicate(4 - length(plots), ggplot() + theme_void(), simplify = FALSE))
  }
  
  combined <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])
  
  #### Save PNG change pattern number
  ggsave(filename = paste0("violin_pattern6_", i, ".png"), plot = combined, width = 10, height = 8)
}


gene= intersect(rownames(cellchat@data), ligand_df$Pattern_3)

plots <- FeaturePlot(seurat_object, features = gene, combine = FALSE)
plots_per_page <- 4
plot_chunks <- split(plots, ceiling(seq_along(plots) / plots_per_page))

# Salva ogni gruppo di grafici in un file PNG separato

for (i in seq_along(plot_chunks)) {
  # Crea il file PNG per ogni pagina
  png(filename = paste0("feature_pattern6_page_", i, ".png"), width = 10, height = 8, units = "in", res = 300)
  
  # Combina i grafici in una griglia (2 colonne)
  gridExtra::grid.arrange(grobs = plot_chunks[[i]], ncol = 2)
  
  # Chiudi il file PNG
  dev.off()
}

gene= intersect(rownames(cellchat@data), ligand_df$Pattern_1)

DotPlot(seurat_object, features = gene,group.by = "pattern")
DoHeatmap(seurat_object, features = gene,group.by = "pattern", size=3)+
  NoLegend()


##################

netAnalysis_signalingRole_network(cellchat, signaling = "SPP1", width = 8, height = 2.5, font.size = 10)
netAnalysis_contribution(cellchat, signaling = "SPP1")
netVisual_heatmap(cellchat, color.heatmap = "Reds", signaling = "SPP1")
netVisual_individual(cellchat, signaling = "SPP1", layout = "circle")

signal_list <- c("SPP1", "MIF", "Adenosine", "FGF", "L1CAM", "EPHB", "EPHA")  # la tua lista di geni

for (signal in signal_list) {
  cat("Plotting:", signal, "\n")
  
  # 1. NETWORK
  png(paste0(signal, "_network.png"), width = 800, height = 250)
  netAnalysis_signalingRole_network(cellchat, signaling = signal, width = 8, height = 2.5, font.size = 10)
  dev.off()
  
  # 2. CONTRIBUTION
  png(paste0(signal, "_contribution.png"), width = 800, height = 800)
  grid::grid.newpage()  # importante per evitare conflitti
  netAnalysis_contribution(cellchat, signaling = signal)
  dev.off()
  
  # 3. HEATMAP
  png(paste0(signal, "_heatmap.png"), width = 800, height = 800)
  grid::grid.newpage()  # anche qui
  netVisual_heatmap(cellchat, color.heatmap = "Reds", signaling = signal)
  dev.off()
}


