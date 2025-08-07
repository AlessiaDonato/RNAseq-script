library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(Orthology.eg.db)

#mouse to human
#tolgo dalla funzione i mouse symbol e mapped
mapIt <- function(mouseids, horg, morg, orth){
  mouseg <- mapIds(morg, mouseids, "ENTREZID", "SYMBOL")
  mapped <- select(orth, mouseg, "Homo_sapiens","Mus_musculus")
  names(mapped) <- c("Mus_egid", "Homo_egid")
  husymb <- select(horg, as.character(mapped[,2]), "SYMBOL","ENTREZID")
  return(data.frame(#Mus_symbol = mouseids,
    #mapped,
    Homo_symbol = husymb[,2]))
}

mapIt("Ado", org.Hs.eg.db, org.Mm.eg.db, 
                             Orthology.eg.db)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(Orthology.eg.db)

#human to mouse
#tolgo dalla funzione i mouse symbol e mapped
mapIt <- function(mouseids, morg, horg, orth){
  homog <- mapIds(horg, mouseids, "ENTREZID", "SYMBOL")
  mapped <- AnnotationDbi::select(orth, homog,"Mus.musculus", "Homo.sapiens")
  names(mapped) <- c("Homo.egid","Mus.egid")
  mosymb <- select(morg, as.character(mapped[,2]), "SYMBOL","ENTREZID")
  return(data.frame(#Mus_symbol = mouseids,
    #mapped,
    mouse.symbol = mosymb[,2]))
}

mapIt("ADO", org.Mm.eg.db, org.Hs.eg.db,
                             Orthology.eg.db)
