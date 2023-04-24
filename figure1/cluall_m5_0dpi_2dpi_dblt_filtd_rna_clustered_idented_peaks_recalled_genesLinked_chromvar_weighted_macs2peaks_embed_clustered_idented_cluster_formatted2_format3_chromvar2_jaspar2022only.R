setwd('/media/shyam/external/multiome_clu_2/multi5')
library(Seurat)
library(Signac)
load("/media/shyam/external/multiome_clu_2/multi5/cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar_weighted_macs2peaks_embed_clustered_idented_cluster_formatted2_format3.RData")
library(JASPAR2022)
library(TFBSTools)
pfm <- getMatrixSet(
x = JASPAR2022,
opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
library(BSgenome.Mmusculus.ensembl.mm39)
DefaultAssay(clu) <- 'peaks'
clu <- AddMotifs(clu, genome = BSgenome.Mmusculus.ensembl.mm39, pfm = pfm)
clu <- RunChromVAR(clu, genome = BSgenome.Mmusculus.ensembl.mm39, new.assay.name = 'chromvar2')
collapse_peaks <- read.csv(file = 'clu_0_2_dpi_dblt_filtd_all_peaks_collapse_clusters_zonationMedian.csv', row.names = 1)
revsc_peaks <- collapse_peaks[collapse_peaks$cluster == 'RevSC1',]
revsc_peaks <- revsc_peaks[revsc_peaks$p_val < 0.005,]
top_peaks_revsc <- revsc_peaks$gene
query <- top_peaks_revsc
DefaultAssay(clu) <- 'peaks'
revsc_motfs2 <- FindMotifs(clu, features = top_peaks_revsc)
write.csv(revsc_motfs2, file = 'clu_0_2_dpi_dblt_filtd_revsc_motifs_new_jaspar2022_list_redo.csv')
DefaultAssay(clu) <- 'RNA'
save(clu, file = 'cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar_weighted_macs2peaks_embed_clustered_idented_cluster_formatted2_format3_chromvar2_jaspar2022only.RData')
savehistory("/media/shyam/external/multiome_clu_2/multi5/cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar_weighted_macs2peaks_embed_clustered_idented_cluster_formatted2_format3_chromvar2_jaspar2022only.R")
