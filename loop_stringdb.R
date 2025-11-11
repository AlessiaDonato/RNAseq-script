library(STRINGdb)
##choose the NCBI taxonomy identifiers of the organism on which you have
#performed the experiment (e.g. 9606 for Human, 10090 for mouse)

string_db <- STRINGdb$new( version="12.0", species=10090, 
                           network_type="full",
                           link_data='combined_only', input_directory="")
#score_threshold:Any interactions with scores below this threshold will not be 
#loaded into the object (the default threshold is 400)
#network_type’: "functional" for the full functional STRING network, or "physical" 
#for the physical subnetwork that links only proteins within the same physical complex.

length(genes_per_module)

for (i in seq_along(genes_per_module)) { 
  cat("Analizzando modulo:", i, "\n")
  
  gene_module <- as.data.frame(genes_per_module[[i]])
  colnames(gene_module) <- "gene"
  
  example1_mapped <- string_db$map(gene_module,"gene",removeUnmappedRows = TRUE)
  write.csv2(example1_mapped, paste0("module_", i, "_mapped.csv"))
  
  hits <- example1_mapped$STRING_id
  enrichment <- string_db$get_enrichment(hits)
  write.csv2(enrichment, paste0("enrichment_module_", i, ".csv"))
  
  # Dividi in chunk di max 2000 se serve
  chunks <- if (length(hits) > 2000) {
    split(hits, ceiling(seq_along(hits) / 2000))
  } else {
    list(hits)
  }
  
  for (j in seq_along(chunks)) {
    cat("  Plot chunk", j,"del modulo", i, "\n")
    png(paste0("string_plot_module_", i, ifelse(length(chunks) > 1, 
                                                paste0("_chunk_", j),""),".png"), 
        width = 3000, height = 2400, res = 300)
    string_db$plot_network(chunks[[j]])
    dev.off()
  }
  
  # Clustering e plot (solo cluster >5)
  clustersList <- string_db$get_clusters(example1_mapped$STRING_id)
  for (a in seq_along(clustersList)) {
    if (length(clustersList[[a]]) > 5) {
      cat("  ➤ Salvando cluster", a, "del modulo", i, "\n")
      file_name <- paste0("module_", i, "_cluster_network_", a, ".png")
      png(filename = file_name, width = 3000, height = 2400, res = 300)
      string_db$plot_network(clustersList[[a]])
      dev.off()
    }
  }
}
