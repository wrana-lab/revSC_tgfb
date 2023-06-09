setwd("~/pCloudDrive/multiome_work/multi5")
library(Seurat)
library(Signac)
load("~/pCloudDrive/multiome_work/multi5/cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar_weighted_macs2peaks_embed_clustered_idented_cluster_formatted2_format3_chromvar2_jaspar2022only.RData")
unique(clu$new_collapse_idents)
Idents(clu) <- clu$new_collapse_idents
clu <- RenameIdents(clu, 'Villus Top' = 'Enterocytes', 'Villus Bottom' = 'Enterocytes', 'CVJ' = 'Enterocytes')
clu$new_collapse_idents2 <- Idents(clu)
DimPlot(clu, reduction = 'wnn.umap')
save(clu, file = 'cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar_weighted_macs2peaks_embed_clustered_idented_cluster_formatted2_format3_chromvar2_jaspar2022only_newcollapse2.RData')
savehistory("~/pCloudDrive/multiome_work/multi5/cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar_weighted_macs2peaks_embed_clustered_idented_cluster_formatted2_format3_chromvar2_jaspar2022only_newcollapse2.R")
