setwd("/media/shyam/external/multiome_clu_2/multi5")
library(Seurat)
library(Signac)
load("/media/shyam/external/multiome_clu_2/multi5/cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar_weighted_macs2peaks_embed.RData")
clu <- FindNeighbors(clu, dims = 2:30, reduction = 'lsi')
clu <- FindClusters(clu, resolution = seq(0.1, 1.5, by = 0.1), graph.name = 'peaks_snn')
library(clustree)
clustree(clu, prefix = 'peaks_snn_res.')
clu <- FindClusters(clu, resolution = 0.5, graph.name = 'peaks_snn')
DimPlot(clu, reduction  = 'umap.peaks')
clu$peaks_clusters <- Idents(clu)
clu$peaks_clusters_idents <- paste('A', clu$peaks_clusters,sep = '' )
Idents(clu) <- clu$peaks_clusters_idents
DimPlot(clu, reduction  = 'umap.peaks')
DefaultAssay(clu) <- 'RNA'
clu <- FindClusters(clu, graph.name = 'wsnn', algorithm = 3, resolution = seq(0.1,2.0,by=0.1))
clustree(clu, prefix = 'wsnn_res.')
clu <- FindClusters(clu, graph.name = 'wsnn', algorithm = 3, resolution = 0.3)
DimPlot(clu, label = T, reduction  = 'wnn.umap')
clu <- FindClusters(clu, graph.name = 'wsnn', algorithm = 3, resolution = 0.5)
DimPlot(clu, label = T, reduction  = 'wnn.umap')
FeaturePlot(clu, features = 'Lgr5', reduction = 'wnn.umap', order = T)
FeaturePlot(clu, features = 'Clu', reduction = 'wnn.umap', order = T)
FeaturePlot(clu, features = 'Atg9b', reduction = 'wnn.umap', order = T)
FeaturePlot(clu, features = 'Defa17', reduction = 'wnn.umap', order = T)
library(ggplot2)
DotPlot(clu, features = c('Alpi', 'Aoc1', 'Ccl25', 'Muc2', 'Agr2', 'Chga', 'Chgb', 'Dclk1', 'Trpm5','Defa17', 'Defa22', 'Lgr5', 'Olfm4', 'Clu', 'F3', 'Atg9b', 'Anxa1', 'Ly6a', 'Mki67','Pdgfra', 'Col1a1', 'Cd3g', 'Cd160', 'Lyz2', 'Csf1r'), col.min = 0) & scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black')) & theme(panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
FeaturePlot(clu, features = c('Alpi', 'Aoc1', 'Ccl25', 'Muc2', 'Agr2', 'Chga', 'Chgb', 'Dclk1', 'Trpm5','Defa17', 'Defa22', 'Lgr5', 'Olfm4', 'Clu', 'F3', 'Atg9b', 'Anxa1', 'Ly6a', 'Mki67','Pdgfra', 'Col1a1', 'Cd3g', 'Cd160', 'Lyz2', 'Csf1r'), reduction = 'umap.rna', order = T, label = T) & scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black'))
clu$wnn_clusters <- Idents(clu)
ids <- read.csv(file = 'clu_0_2_dpi_dblt_filtd_wnn_clusters_idents.csv', row.names = 1)
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
DimPlot(clu, label = T, reduction = 'wnn.umap')
FeaturePlot(clu, features = c('Alpi', 'Aoc1', 'Ccl25', 'Muc2', 'Agr2', 'Chga', 'Chgb', 'Dclk1', 'Trpm5','Defa17', 'Defa22', 'Lgr5', 'Olfm4', 'Clu', 'F3', 'Atg9b', 'Anxa1', 'Ly6a', 'Mki67','Pdgfra', 'Col1a1', 'Cd3g', 'Cd160', 'Lyz2', 'Csf1r'), reduction = 'wnn.umap', order = T, label = T) & scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black'))
Idents(clu) <- clu$wnn_clusters
ids <- read.csv(file = 'clu_0_2_dpi_dblt_filtd_wnn_clusters_idents.csv', row.names = 1)
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
DimPlot(clu, label = T, reduction = 'wnn.umap')
clu$wnn_clusters_idents <- Idents(clu)
FeaturePlot(clu, features = c('Dclk1', 'Trpm5'), reduction = 'wnn.umap', order = T)
FeaturePlot(clu, features = c('Dclk1'), reduction = 'wnn.umap', order = T)
p <- FeaturePlot(clu, features = c('Dclk1'), reduction = 'wnn.umap', order = T)
sel <- CellSelector(p)
clu <- SetIdent(clu, cells = sel, value = 'TC')
clu$wnn_clusters_idents <- Idents(clu)
Idents(clu) <- clu$wnn_clusters_idents
DimPlot(clu, reduction = 'wnn.umap')
DimPlot(clu, reduction = 'wnn.umap', label = T)
save(clu, file = 'cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar_weighted_macs2peaks_embed_clustered_idented.RData')
savehistory("/media/shyam/external/multiome_clu_2/multi5/cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar_weighted_macs2peaks_embed_clustered_idented.R")
