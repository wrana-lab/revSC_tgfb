setwd("/media/shyam/external/multiome_clu_2/multi5")
library(Seurat)
library(Signac)
load("/media/shyam/external/multiome_clu_2/multi5/cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar_weighted_macs2peaks_embed_clustered_idented_cluster_formatted2_format3_chromvar2_jaspar2022only_newcollapse2_rezonationed_mean.RData")
collapse_peaks <- read.csv(file = 'clu_0_2_dpi_dblt_filtd_peaks_collapse_zonationmean_zoneswise.csv', row.names = 1)
DefaultAssay(clu) <- 'peaks'
annos <- Annotation(clu)
annos <- as.data.frame(annos)
total_peaks <- rownames(clu)
test <- StringToGRanges(total_peaks)
library(BSgenome.Mmusculus.ensembl.mm39)
mouse_views2 <- BSgenomeViews(BSgenome.Mmusculus.ensembl.mm39, test)
library(JASPAR2022)
library(TFBSTools)
pfm <- getMatrixSet(
x = JASPAR2022,
opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
tfs <- as.data.frame(t(as.data.frame(clu@assays$peaks@motifs@motif.names)))
names(pfm) <- tfs[names(pfm),]
pfm_query <- pfm
for(i in names(pfm_query)){
if(i %in% c('TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'SMAD2')){
next
} else {
pfm_query[[i]] <- NULL
}
}
library(motifmatchr)
peaks_mat <- clu@assays$peaks@counts
revsc_peaks <- collapse_peaks[collapse_peaks$cluster == 'RevSC1',]
revsc_peaks <- revsc_peaks[revsc_peaks$p_val < 0.005,]
revsc_peaks <- collapse_peaks[collapse_peaks$cluster == 'RevSC1-W13',]
revsc_peaks <- revsc_peaks[revsc_peaks$p_val < 0.005,]
top_peaks_revsc <- revsc_peaks$gene
query <- top_peaks_revsc
df <- StringToGRanges(query)
mouse_motifs2 <- matchMotifs(pwms = pfm_query, subject = df, out = c('matches'), genome = BSgenome.Mmusculus.ensembl.mm39, bg = colSums(letterFrequency(mouse_views2, letters = c('A', 'C', 'G', 'T'), as.prob = F)) / sum(colSums(letterFrequency(mouse_views2, letters = c('A', 'C', 'G', 'T'), as.prob = F))))
###########
matches_per_peak <- as.matrix(motifMatches(mouse_motifs2))
matches_per_peak <- as.data.frame(matches_per_peak)
matches_per_peak$TEAD <- ifelse(matches_per_peak$TEAD3 == T | matches_per_peak$TEAD2 == T | matches_per_peak$TEAD1 == T | matches_per_peak$TEAD4 == T, T, F)
matches_per_peak$smad <- ifelse(matches_per_peak$SMAD2 == T, T, F)
View(matches_per_peak)
matches_per_peak$TEAD_only <- ifelse(matches_per_peak$TEAD == T & matches_per_peak$smad == F, T, F)
matches_per_peak$smad_only <- ifelse(matches_per_peak$TEAD == F & matches_per_peak$smad == T, T, F)
matches_per_peak$smad_tead_both <- ifelse(matches_per_peak$TEAD == T & matches_per_peak$smad == T, T, F)
matches_per_peak_revsc <- matches_per_peak[,c('TEAD_only', 'smad_only', 'smad_tead_both')]
matches_per_peak_revsc$peaks <- top_peaks_revsc
View(matches_per_peak_revsc)
df_total <- rownames(peaks_mat)
df_total <- StringToGRanges(df_total)
mouse_motifstotal <- matchMotifs(pwms = pfm_query, subject = df_total, out = c('matches'), genome = BSgenome.Mmusculus.ensembl.mm39, bg = colSums(letterFrequency(mouse_views2, letters = c('A', 'C', 'G', 'T'), as.prob = F)) / sum(colSums(letterFrequency(mouse_views2, letters = c('A', 'C', 'G', 'T'), as.prob = F))))
#matches per every peak in global peakset
matches_per_peak_total <- as.matrix(motifMatches(mouse_motifstotal))
matches_per_peak_total <- as.data.frame(matches_per_peak_total)
matches_per_peak_total$TEAD <- ifelse(matches_per_peak_total$TEAD3 == T | matches_per_peak_total$TEAD2 == T | matches_per_peak_total$TEAD1 == T | matches_per_peak_total$TEAD4 == T, T, F)
matches_per_peak_total$smad <- ifelse( matches_per_peak_total$SMAD2 == T, T, F)
matches_per_peak_total$TEAD_only <- ifelse(matches_per_peak_total$TEAD == T & matches_per_peak_total$smad == F, T, F)
matches_per_peak_total$smad_only <- ifelse(matches_per_peak_total$TEAD == F & matches_per_peak_total$smad == T, T, F)
matches_per_peak_total$smad_tead_both <- ifelse(matches_per_peak_total$TEAD == T & matches_per_peak_total$smad == T, T, F)
matches_per_peak_global <- matches_per_peak_total[,c('TEAD_only', 'smad_only', 'smad_tead_both')]
matches_per_peak_global$peaks <- rownames(peaks_mat)
View(matches_per_peak_global)
write.csv(matches_per_peak_global, file = 'global_peakset_smad_tead_both_clu_0_2_dpi_dblt_filtd.csv')
write.csv(matches_per_peak_revsc, file = 'revsc_selective_peakset_smad_tead_both_clu_0_2_dpi_dblt_filtd.csv')
savehistory("/media/shyam/external/multiome_clu_2/multi5/global_peakset_smad_tead_both_clu_0_2_dpi_dblt_filtd.R")
