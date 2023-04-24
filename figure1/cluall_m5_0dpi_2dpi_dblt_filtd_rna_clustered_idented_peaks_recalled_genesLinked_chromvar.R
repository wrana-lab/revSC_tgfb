setwd("/media/shyam/external/multiome_clu_2/multi5")
library(Seurat)
library(Signac)
load("/media/shyam/external/multiome_clu_2/multi5/cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked.RData")
library(BSgenome.Mmusculus.ensembl.mm39)
library(JASPAR2020)
library(TFBSTools)
DefaultAssay(clu) <- 'peaks'
pfm <- getMatrixSet(
x = JASPAR2020,
opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
clu <- AddMotifs(clu, genome = BSgenome.Mmusculus.ensembl.mm39, pfm = pfm)
clu <- RunChromVAR(clu, genome = BSgenome.Mmusculus.ensembl.mm39)
save(clu, file = 'cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar.RData')
savehistory("/media/shyam/external/multiome_clu_2/multi5/cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar.R")
