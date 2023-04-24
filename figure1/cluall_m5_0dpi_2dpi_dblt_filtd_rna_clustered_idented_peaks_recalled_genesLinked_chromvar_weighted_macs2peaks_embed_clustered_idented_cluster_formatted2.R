setwd("/media/shyam/external/multiome_clu_2/multi5")
library(Seurat)
library(Signac)
load("/media/shyam/external/multiome_clu_2/multi5/cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar_weighted_macs2peaks_embed_clustered_idented_cluster_formatted.RData")
Idents(clu) <- clu$new_wnn_clusters
motiffss <- as.data.frame(t(as.data.frame(clu@assays$peaks@motifs@motif.names)))
DefaultAssay(clu) <- 'chromvar'
rownames(clu@assays$chromvar@data) <- motiffss$V1
wnn_clusters_motifs <- FindAllMarkers(clu, only.pos = T, mean.fxn = rowMeans, fc.name= 'avg_diff')
library(dplyr)
wnn_clusters_motifs %>% group_by(cluster) %>% top_n(n = 25, wt = avg_diff) -> top25
write.csv(top25, file = 'clu_0_2dpi_dblt_filtd_wnnclusters_top25_motifs.csv')
Idents(clu) <- clu$peaks_clusters_idents
peaks_clusters_motifs <- FindAllMarkers(clu, only.pos = T, mean.fxn = rowMeans, fc.name= 'avg_diff')
peaks_clusters_motifs %>% group_by(cluster) %>% top_n(n = 25, wt = avg_diff) -> top25_peaks
Idents(clu) <- clu$rna_clusters_idents
for(i in c('EC', 'GC', 'EE', 'PC', 'TC','ICC' ,'Immune', 'Fibro')){
clu <- SetIdent(clu, cells = WhichCells(clu, expression = rna_clusters_idents %in% grep(paste('^', i, '.*', sep = ""), unique(clu$rna_clusters_idents), value = T)), value = i)
}
clu$new_rna_clusters <- paste(Idents(clu), '-','R', clu$rna_clusters, sep = "")
DimPlot(clu, reduction = 'wnn.umap', label = T)
Idents(clu) <- clu$new_rna_clusters
DimPlot(clu, reduction = 'wnn.umap', label = T)
rna_clusters_motifs <- FindAllMarkers(clu, only.pos = T, mean.fxn = rowMeans, fc.name= 'avg_diff')
rna_clusters_motifs %>% group_by(cluster) %>% top_n(n = 25, wt = avg_diff) -> top25_rna
View(wnn_clusters_motifs)
write.csv(top25_peaks, file = 'clu_0_2dpi_dblt_filtd_peaksclusters_top25_motifs.csv')
write.csv(top25_rna, file = 'clu_0_2dpi_dblt_filtd_rnaclusters_top25_motifs.csv')
chrom_mat <- as.matrix(clu@assays$chromvar@data)
library(RcppCNPy)
npySave(filename = 'clu_0_2_dpi_dblt_filtd_chomvar_mat.npy', chrom_mat)
motif_names <- as.data.frame(rownames(chrom_mat))
View(motif_names)
write.csv(motif_names, file = 'clu_0_2_dpi_dblt_filtd_chromv_motif_names.csv')
gc()
peaksmat <- clu@assays$peaks@data
library(reticulate)
np <- import('numpy')
np$save(file = 'clu_0_2_dpi_dblt_filtd_peaks_mat.npy', peaksmat)
np$save(peaksmat, file = 'clu_0_2_dpi_dblt_filtd_peaks_mat.npy')
np$save('clu_0_2_dpi_dblt_filtd_peaks_mat.npy', peaksmat)
peaksnames <- as.data.frame(rownames(peaksmat))
write.csv(peaksnames, file = 'clu_0_2_dpi_dblt_filtd_peaksnames.csv')
meta <- clu@meta.data
write.csv(meta, file= 'clu_0_2_dpi_dblt_filtd_meta.csv')
DotPlot(clu, features = c('FOS', 'FOS::JUN', 'JUN(var.2)', 'SMAD3', 'Smad4', 'Smad2::Smad3', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4'))
Idents(clu) <- clu$new_wnn_clusters
Idents(clu) <- clu$collapse_idents
DotPlot(clu, features = c('FOS', 'FOS::JUN', 'JUN(var.2)', 'SMAD3', 'Smad4', 'Smad2::Smad3', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4'))
DotPlot(clu, features = c('FOS', 'FOS::JUN', 'JUN(var.2)', 'SMAD3', 'Smad4', 'Smad2::Smad3', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4'), idents = levels(clu)[!levels(clu) %in% c(grep('Immune', levels(clu), value = T), grep('Fibro', levels(clu), value = T))])
levels(clu) <- c(grep('^RevSC.*', levels(clu), value = T) ,grep('^FCC.*', levels(clu), value = T),grep('^Crypt.*', levels(clu), value = T), grep('^Paneth.*', levels(clu), value = T),grep('^Tuft.*', levels(clu), value = T),grep('^Enteroendo.*', levels(clu), value = T), grep('^Goblet.*', levels(clu), value = T),grep('^Enterocytes.*', levels(clu), value = T), grep('^31.*', levels(clu), value = T),grep('^Fibro.*', levels(clu), value = T),grep('^Immune.*', levels(clu), value = T))
DotPlot(clu, features = c('FOS', 'FOS::JUN', 'JUN(var.2)', 'SMAD3', 'Smad4', 'Smad2::Smad3', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4'), idents = levels(clu)[!levels(clu) %in% c(grep('Immune', levels(clu), value = T), grep('Fibro', levels(clu), value = T))])
levels(clu) <- c(grep('^RevSC.*', levels(clu), value = T) ,grep('^FCC.*', levels(clu), value = T),grep('^Crypt.*', levels(clu), value = T),grep('^Enterocytes.*', levels(clu), value = T), grep('^Paneth.*', levels(clu), value = T),grep('^Tuft.*', levels(clu), value = T),grep('^Enteroendo.*', levels(clu), value = T), grep('^Goblet.*', levels(clu), value = T),grep('^31.*', levels(clu), value = T),grep('^Fibro.*', levels(clu), value = T),grep('^Immune.*', levels(clu), value = T))
DotPlot(clu, features = c('FOS', 'FOS::JUN', 'JUN(var.2)', 'SMAD3', 'Smad4', 'Smad2::Smad3', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4'), idents = levels(clu)[!levels(clu) %in% c(grep('Immune', levels(clu), value = T), grep('Fibro', levels(clu), value = T))])
DotPlot(clu, features = c('FOS', 'FOS::JUN', 'JUN(var.2)', 'SMAD3', 'Smad4', 'Smad2::Smad3', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4'), idents = levels(clu)[!levels(clu) %in% c(grep('Immune', levels(clu), value = T), grep('Fibro', levels(clu), value = T))]) & scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black')) & theme(panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
library(ggplot2)
DotPlot(clu, features = c('FOS', 'FOS::JUN', 'JUN(var.2)', 'SMAD3', 'Smad4', 'Smad2::Smad3', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4'), idents = levels(clu)[!levels(clu) %in% c(grep('Immune', levels(clu), value = T), grep('Fibro', levels(clu), value = T))]) & scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black')) & theme(panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(clu, features = c('FOS', 'FOS::JUN', 'JUN(var.2)', 'SMAD3', 'Smad4', 'Smad2::Smad3', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4'), idents = levels(clu)[!levels(clu) %in% c(grep('Immune', levels(clu), value = T), grep('Fibro', levels(clu), value = T))]) & scale_color_gradientn(colors = c('#2166ac', '#b2182b')) & theme(panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(clu, features = c('FOS', 'FOS::JUN', 'JUN(var.2)', 'SMAD3', 'Smad4', 'Smad2::Smad3', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4'), idents = levels(clu)[!levels(clu) %in% c(grep('Immune', levels(clu), value = T), grep('Fibro', levels(clu), value = T))]) & scale_color_gradientn(colors = c('#2166ac','#ffffff' ,'#b2182b')) & theme(panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(clu, features = c('FOS', 'FOS::JUN', 'JUN(var.2)', 'SMAD3', 'Smad4', 'Smad2::Smad3', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4'), idents = levels(clu)[!levels(clu) %in% c(grep('Immune', levels(clu), value = T), grep('Fibro', levels(clu), value = T))]) & scale_color_gradientn(colors = c('#2166ac','#ffffff' ,'#b2182b'), midpoint = 0) & theme(panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(clu, features = c('FOS', 'FOS::JUN', 'JUN(var.2)', 'SMAD3', 'Smad4', 'Smad2::Smad3', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4'), idents = levels(clu)[!levels(clu) %in% c(grep('Immune', levels(clu), value = T), grep('Fibro', levels(clu), value = T))]) & scale_color_gradient2(colors = c('#2166ac','#ffffff' ,'#b2182b'), midpoint = 0) & theme(panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(clu, features = c('FOS', 'FOS::JUN', 'JUN(var.2)', 'SMAD3', 'Smad4', 'Smad2::Smad3', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4'), idents = levels(clu)[!levels(clu) %in% c(grep('Immune', levels(clu), value = T), grep('Fibro', levels(clu), value = T))]) & scale_color_gradient2(low = '#2166ac', mid = '#ffffff', high ='#b2182b'), midpoint = 0) & theme(panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(clu, features = c('FOS', 'FOS::JUN', 'JUN(var.2)', 'SMAD3', 'Smad4', 'Smad2::Smad3', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4'), idents = levels(clu)[!levels(clu) %in% c(grep('Immune', levels(clu), value = T), grep('Fibro', levels(clu), value = T))]) & scale_color_gradient2(low = '#2166ac', mid = '#ffffff', high ='#b2182b', midpoint = 0) & theme(panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DefaultAssay(clu)
DefaultAssay(clu) <- 'RNA'
smads <- c('Cdkn2b', 'Cdkn1a', 'Fn1', 'Col1a1', 'Col3a1', 'Col5a1', 'Col6a1',
'Dcn', 'Bgn', 'Acan', 'Spock1', 'Spock2', 'Spock3', 'Hspg2',
'Tgfbr3', 'Agrn', 'Ambp', 'Ncan', 'Vcan', 'Bcan', 'Fmod', 'Lum',
'Tnc', 'Tnr', 'Tnxb', 'Sparc', 'Spp1', 'Thbs1', 'Timp1', 'Timp2',
'Timp3', 'Timp4', 'Vtn', 'Serpine1', 'Itga1', 'Itga2', 'Itga3',
'Itga4', 'Itga5', 'Itga6', 'Itga7', 'Itga10', 'Itga11', 'Itgav',
'Itga2b', 'Itgb1', 'Itgb2', 'Itgb3', 'Itgb4', 'Itgb5', 'Itgb6',
'Itgb8', 'Thbs2', 'Thbs3', 'Thbs4', 'Eln', 'Mmp1b', 'Mmp1a',
'Mmp8', 'Mmp13', 'Mmp3', 'Itgb2l')
clu$smads <- apply(clu@assays$RNA@data[smads[smads %in% rownames(clu@assays$RNA@data)],], 2, mean)
teads <- c('Crim1', 'Ccn1', 'Ccn2', 'Amotl2', 'Ankrd1', 'Igfbp3', 'F3',
'Fjx1', 'Nuak2', 'Lats2', 'Gm49361', 'Gadd45a', 'Tgfb2', 'Ptpn14',
'Nt5e', 'Foxf2', 'Axl', 'Dock5', 'Asap1', 'Rbms3', 'Myof',
'Arhgef17', 'Ccdc80')
clu$teads <- apply(clu@assays$RNA@data[teads[teads %in% rownames(clu@assays$RNA@data)],], 2, mean)
DotPlot(clu, features = c('smads', 'teads'), col.min = 0, idents= c('Enterocytes', 'Goblet Cells', 'Enteroendocrine Cells', 'Paneth Cells', 'Crypt Cells', 'RevSC1')) & ggplot2::scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black')) &  labs(x = 'Pathway Expression') & ggplot2::theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15))
clu$smads <- apply(clu@assays$RNA@data[smads[smads %in% rownames(clu@assays$RNA@data)],], 2, median)
teads <- c('Crim1', 'Ccn1', 'Ccn2', 'Amotl2', 'Ankrd1', 'Igfbp3', 'F3',
'Fjx1', 'Nuak2', 'Lats2', 'Gm49361', 'Gadd45a', 'Tgfb2', 'Ptpn14',
'Nt5e', 'Foxf2', 'Axl', 'Dock5', 'Asap1', 'Rbms3', 'Myof',
'Arhgef17', 'Ccdc80')
clu$teads <- apply(clu@assays$RNA@data[teads[teads %in% rownames(clu@assays$RNA@data)],], 2, median)
DotPlot(clu, features = c('smads', 'teads'), col.min = 0, idents= c('Enterocytes', 'Goblet Cells', 'Enteroendocrine Cells', 'Paneth Cells', 'Crypt Cells', 'RevSC1')) & ggplot2::scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black')) &  labs(x = 'Pathway Expression') & ggplot2::theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15))
smads <- c('Cdkn2b', 'Cdkn1a', 'Fn1', 'Col1a1', 'Col3a1', 'Col5a1', 'Col6a1',
'Dcn', 'Bgn', 'Acan', 'Spock1', 'Spock2', 'Spock3', 'Hspg2',
'Tgfbr3', 'Agrn', 'Ambp', 'Ncan', 'Vcan', 'Bcan', 'Fmod', 'Lum',
'Tnc', 'Tnr', 'Tnxb', 'Sparc', 'Spp1', 'Thbs1', 'Timp1', 'Timp2',
'Timp3', 'Timp4', 'Vtn', 'Serpine1', 'Itga1', 'Itga2', 'Itga3',
'Itga4', 'Itga5', 'Itga6', 'Itga7', 'Itga10', 'Itga11', 'Itgav',
'Itga2b', 'Itgb1', 'Itgb2', 'Itgb3', 'Itgb4', 'Itgb5', 'Itgb6',
'Itgb8', 'Thbs2', 'Thbs3', 'Thbs4', 'Eln', 'Mmp1b', 'Mmp1a',
'Mmp8', 'Mmp13', 'Mmp3', 'Itgb2l')
clu$smads <- apply(clu@assays$RNA@data[smads[smads %in% rownames(clu@assays$RNA@data)],], 2, mean)
teads <- c('Crim1', 'Ccn1', 'Ccn2', 'Amotl2', 'Ankrd1', 'Igfbp3', 'F3',
'Fjx1', 'Nuak2', 'Lats2', 'Gm49361', 'Gadd45a', 'Tgfb2', 'Ptpn14',
'Nt5e', 'Foxf2', 'Axl', 'Dock5', 'Asap1', 'Rbms3', 'Myof',
'Arhgef17', 'Ccdc80')
clu$teads <- apply(clu@assays$RNA@data[teads[teads %in% rownames(clu@assays$RNA@data)],], 2, mean)
DotPlot(clu, features = c('smads', 'teads'), col.min = 0, idents= c('Enterocytes', 'Goblet Cells', 'Enteroendocrine Cells', 'Paneth Cells', 'Crypt Cells', 'RevSC1')) & ggplot2::scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black')) &  labs(x = 'Pathway Expression') & ggplot2::theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15))
save(clu, file = 'cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar_weighted_macs2peaks_embed_clustered_idented_cluster_formatted2.RData')
savehistory("/media/shyam/external/multiome_clu_2/multi5/cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar_weighted_macs2peaks_embed_clustered_idented_cluster_formatted2.R")
