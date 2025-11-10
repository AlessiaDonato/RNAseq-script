Convert("immune_noID_seurat.h5ad", dest = "h5seurat", overwrite = F)
immune<- LoadH5Seurat("immune_noID_seurat.h5seurat", assays="RNA", meta.data = FALSE, misc = FALSE)
obs <- h5read("immune_noID_seurat.h5ad", "/obs")
meta <- data.frame(lapply(names(obs), function(x) { 
  if (length(obs[[x]])==2) 
    obs[[x]][['categories']][ifelse(obs[[x]][['codes']] >= 0, obs[[x]][['codes']] + 1, NA)]
  else 
    as.numeric(obs[[x]])
}
), row.names=Cells(immune))
colnames(meta) <- names(obs)

immune <- AddMetaData(immune,meta)