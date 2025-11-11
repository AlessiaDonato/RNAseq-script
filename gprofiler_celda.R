library(zellkonverter)
library(SingleCellExperiment)
library(celda)
library(scater)
library(clusterProfiler)
library(gprofiler2)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(ggplot2)

sce<- readRDS("celda_5counts_50cell.rds")
png("celda_heatmap_hd.png", width = 4000, height = 4000, res = 300)
plot(celdaHeatmap(sce = sce, nfeatures = 10, showNamesFeature=T))
dev.off()

## Extract list of genes assigned to each module
modules <- celdaModules(sce, altExpName ="featureSubset")
names(modules) <- rownames(altExp(sce, "featureSubset"))
genes_per_module <- split(names(modules), modules)

#compute enrichment analisys and plot 
#gost performs statistical enrichment analysis to find over-representation 
#of functions from Gene Ontology, biological pathways like KEGG and Reactome, 
#human disease annotations, etc. 
#This is done with the hypergeometric test 
#followed by correction for multiple testing

for (mod in  1: length( genes_per_module)) {
  results <- gost(query = genes_per_module[mod],
                  organism = "mmusculus", correction_method = "fdr", highlight=T
                  #sources = c("KEGG", "REAC", "GO:BP") 
  )
  df<- results[["result"]]
  df$parents <- sapply(df$parents, toString)  
  df<-df[order(df$p_value),]
  write.csv(df, paste0("modules_", mod, ".csv"), row.names = FALSE)
  
  publish_gosttable(df,highlight_terms = df[1:20,],
                    use_colors = TRUE, 
                    show_columns = c("source", "term_name", "term_size", 
                                     "intersection_size"),
                    filename = paste0("table_modules_", mod, ".png"))
  publish_gosttable(df,highlight_terms = df[df$highlighted==TRUE,],
                    use_colors = TRUE, 
                    show_columns = c("source", "term_name", "term_size", 
                                     "intersection_size"),
                    filename = paste0("table_modules_Highlight", mod, ".png"))
  # heatmap
  png(paste0("plot_module_", mod, ".png"), width = 3000, height =2000, res =300)
  print(gostplot(results, capped = F, interactive = F))
  dev.off()
  gostplot_res<-gostplot(results, capped = F, interactive = T)
  htmlwidgets::saveWidget(gostplot_res, 
                          paste0("plot_module_", mod, "_interattivo.html"),
                          selfcontained = TRUE)
}



