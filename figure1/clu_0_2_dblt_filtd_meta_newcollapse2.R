setwd("~/pCloudDrive/multiome_work/multi5")
library(Seurat)
library(Signac)
load("~/pCloudDrive/multiome_work/multi5/cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar_weighted_macs2peaks_embed_clustered_idented_cluster_formatted2_format3_chromvar2_jaspar2022only_newcollapse2.RData")
Idents(clu) <- clu$new_collapse_idents2
levels(clu)
DefaultAssay(clu) <- 'peaks'
DimPlot(clu, reduction = 'wnn.umap', label = T)
meta <- clu@meta.data
write.csv(meta, file = 'clu_0_2_dblt_filtd_meta_newcollapse2.csv')
savehistory("~/pCloudDrive/multiome_work/multi5/clu_0_2_dblt_filtd_meta_newcollapse2.R")
