######vennplot

pattern_1_4=read.csv("Cellchat/pattern1_vs_pattern4.csv", row.names = 1)
pattern_4_6=read.csv("Cellchat/pattern4_vs_pattern6.csv", row.names = 1)
pattern_6_3=read.csv("Cellchat/pattern6_vs_pattern3.csv", row.names = 1)

KO_ICvs_DCT2<-read.csv("ko_ICtransfvsDCT2.csv", row.names = 1)
KO_ICvs_DCT2_d<- KO_ICvs_DCT2[KO_ICvs_DCT2$KO_IC_trans_pvals_adj < 0.05 &
                                 KO_ICvs_DCT2$KO_IC_trans_logfoldchanges < -1, ]
ranks <-as.numeric(KO_ICvs_DCT2_d$KO_IC_trans_logfoldchanges)
names(ranks) <- KO_ICvs_DCT2_d$KO_IC_trans_names
ranks<- sort(ranks, decreasing = F)

new0vs1=read.csv("leidencom_0vs1_bis.csv", row.names = 1)

intersect_tot<- intersect(intersect(KO_ICvsWT$KO_IC_trans_names
                                    [KO_ICvsWT$KO_IC_trans_pvals_adj<0.05], 
                                    KO_ICvs_DCT2$KO_IC_trans_names
                                    [KO_ICvs_DCT2$KO_IC_trans_pvals_adj<0.05]),
                          new0vs1$X0_names[new0vs1$X0_logfoldchanges<0.05])

library(VennDiagram)

venn.plot <- venn.diagram(
  x = list(lista1 = pattern_1_4$names[pattern_1_4$pvals_adj<0.05], 
           lista2 = pattern_4_6$names[pattern_4_6$pvals_adj<0.05], 
           lista3 = pattern_6_3$names[pattern_6_3$pvals_adj<0.05]),
  category.names = c("pattern_1vs4", "pattern_4vs6", 
                     "pattern_6vs3"),
  filename = NULL,  # Imposta NULL per visualizzare il diagramma direttamente
  output = TRUE,
  fill = c("lightblue", "lightgreen", "lightcoral"),  # Colori dei set
  cat.col = c("blue", "darkgreen", "red"),  # Colore delle etichette
  cat.cex = 1,  # Dimensione del testo delle etichette
  cex = 1,  # Dimensione del testo nei cerchi
  # Impostazioni di sfondo
  plot.background = "white",  # Colore di sfondo (bianco)
  transparent = FALSE  # Impedisce la trasparenza
)

grid.draw(venn.plot)

library(enrichR)
enrichR_results <- enrichr(intersect_tot, databases = c("WikiPathways_2024_Mouse", 
                                                        "KEGG_2019_Mouse",
                                                        "HDSigDB_Mouse_2021"))

plotEnrich(enrichR_results[["WikiPathways_2024_Mouse"]], showTerms = 20, numChar = 40, 
           y = "Count", orderBy = "P.value", title="WikiPathways_2024_Mouse")

# Visualizza i risultati
print(results)
###fgsea
Kegg <- msigdbr("mouse", category="C2", subcategory = "CP:KEGG")
pathways_kegg <- split(Kegg$gene_symbol, Kegg$gs_name)

wiki <- msigdbr("mouse", category="C2", subcategory = "CP:WIKIPATHWAYS")
pathways_wiki <- split(wiki$gene_symbol, wiki$gs_name)

GO_BP <- msigdbr("mouse", category="C5", subcategory = "GO:BP")
pathways_bp <- split(GO_BP$gene_symbol, GO_BP$gs_name)

KO_ICvsWT=read.csv("ko_ICtransfvsICwt.csv", row.names = 1)
dim(KO_ICvsWT)
#[1] 27836    10
KO_ICvsWT_d<- KO_ICvsWT[KO_ICvsWT$KO_IC_trans_pvals_adj< 0.05 &
                           KO_ICvsWT$KO_IC_trans_logfoldchanges < -1.5,1:5]
dim(KO_ICvsWT_d)
#[1] 177   5
ranks <-as.numeric(KO_ICvsWT_d$KO_IC_trans_logfoldchanges)
names(ranks) <- KO_ICvsWT_d$KO_IC_trans_names
ranks<- sort(ranks, decreasing = F)
pathway_list<-list(pathways_bp= pathways_bp, pathways_kegg=pathways_kegg, 
                   pathways_wiki=pathways_wiki)

for (i in 1:length(pathway_list)) {
  # Estrai il pathway corrente
  pathway <- pathway_list[[3]]
  fgsea_results <- fgsea(pathways = pathway, stats = ranks)
  # Ordina i risultati per p-value
  filtered_results <- fgsea_results %>%
    filter(padj < 0.05) %>%
    arrange(padj)
  # Salva i risultati nella lista
  fwrite(filtered_results, paste0("./KO_ICvsDCT2/", names(pathway_list)[i], 
                                  "_fgsea.csv"))
  #results_list[[names(pathway_list)[i]]] <- filtered_results
}


