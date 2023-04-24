setwd("/media/shyam/external/multiome_clu_2/multi5")
library(Seurat)
library(Signac)
load("/media/shyam/external/multiome_clu_2/multi5/cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered.RData")
library(clustree)
clustree(clu, prefix = 'RNA_snn.')
clustree(clu, prefix = 'RNA_snn_res.')
clu <- FindClusters(clu, resolution=seq(0.1,2.0, by = 0.1), graph.name = 'RNA_snn')
clustree(clu, prefix = 'RNA_snn_res.')
clu <- FindClusters(clu, resolution = 1.3)
clu <- FindClusters(clu, resolution = 1.3,  graph.name = 'RNA_snn')
DimPlot(clu, reduction = 'umap.rna', label = T)
clu$rna_clusters <- Idents(clu)
library(ggplot2)
DotPlot(clu, features = c('Alpi', 'Aoc1', 'Ccl25', 'Muc2', 'Agr2', 'Chga', 'Chgb', 'Dclk1', 'Trpm5','Defa17', 'Defa22', 'Lgr5', 'Olfm4', 'Clu', 'F3', 'Atg9b', 'Anxa1', 'Ly6a', 'Mki67','Pdgfra', 'Col1a1', 'Cd3g', 'Cd160', 'Lyz2', 'Csf1r'), col.min = 0) & scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black')) & theme(panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
FeaturePlot(clu, features = c('Alpi', 'Aoc1', 'Ccl25', 'Muc2', 'Agr2', 'Chga', 'Chgb', 'Dclk1', 'Trpm5','Defa17', 'Defa22', 'Lgr5', 'Olfm4', 'Clu', 'F3', 'Atg9b', 'Anxa1', 'Ly6a', 'Mki67','Pdgfra', 'Col1a1', 'Cd3g', 'Cd160', 'Lyz2', 'Csf1r'), reduction = 'umap.rna', order = T, label = T) & scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black'))
FeaturePlot(clu, features = c('Dclk1', 'Trpm5'), reduction = 'umap.rna', order = T, label = T) & scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black'))
FeaturePlot(clu, features = c('Dclk1'), reduction = 'umap.rna', order = T, label = T) & scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black'))
p <- FeaturePlot(clu, features = c('Dclk1'), reduction = 'umap.rna', order = T, label = T) & scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black'))
sel <- CellSelector(p)
sel
ids <- read.csv(file = 'clu_0_2_dpi_dblt_filtd_idents.csv', row.names = 1)
ogclus <- levels(clu)
clusterer <- function(ogclus, ids)  {
require(stringr)
tes = rep(NA, length(ogclus))
for(i in rownames(ids)){
x = str_split(ids[i,], ',')
for(ii in 1:length(x[[1]])){
tes[match(x[[1]][ii], ogclus)] = paste0(i,ii, sep='')
}
}
xx = ogclus[is.na(tes)]
for(i in 1:length(xx)){
tes[match(xx[i], ogclus)] = xx[i]
}
return(tes)
}
newids <- clusterer(ogclus, ids)
names(newids) <- levels(clu)
clu <- RenameIdents(clu, newids)
DimPlot(clu, label = T, reduction = 'umap.rna')
Idents(clu) <- clu$rna_clusters
twenty5 <- FindMarkers(clu, ident.1 = 25, only.pos = T)
View(twenty5)
ids <- read.csv(file = 'clu_0_2_dpi_dblt_filtd_idents.csv', row.names = 1)
ogclus <- levels(clu)
clusterer <- function(ogclus, ids)  {
require(stringr)
tes = rep(NA, length(ogclus))
for(i in rownames(ids)){
x = str_split(ids[i,], ',')
for(ii in 1:length(x[[1]])){
tes[match(x[[1]][ii], ogclus)] = paste0(i,ii, sep='')
}
}
xx = ogclus[is.na(tes)]
for(i in 1:length(xx)){
tes[match(xx[i], ogclus)] = xx[i]
}
return(tes)
}
newids <- clusterer(ogclus, ids)
names(newids) <- levels(clu)
clu <- RenameIdents(clu, newids)
DimPlot(clu, label = T, reduction = 'umap.rna')
clu <- SetIdent(clu, cells = sel, value = 'TC')
DimPlot(clu, label = T, reduction = 'umap.rna')
clu$rna_clusters_idents <- Idents(clu)
save(clu, file = 'cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented.RData')
savehistory("/media/shyam/external/multiome_clu_2/multi5/cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented.R")
