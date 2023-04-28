setwd("/media/shyam/external/multiome_clu_2/multi5")
library(Seurat)
library(Signac)
load("cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar_weighted_macs2peaks_embed_clustered_idented_cluster_formatted2_format3_chromvar2_jaspar2022only_newcollapse2_rezonationed_mean.RData")
Idents(clu) <- clu$new_collapse_idents
levels(clu) <- c(grep('^RevSC', levels(clu), value = T) ,grep('^FCC', levels(clu), value = T),grep('^Crypt Cells', levels(clu), value = T), grep('^Paneth Cells', levels(clu), value = T),grep('^CVJ', levels(clu), value = T),grep('^Villus Bottom', levels(clu), value = T), grep('^Villus Middle', levels(clu), value = T), grep('^Villus Top', levels(clu), value = T), grep('^Goblet Cells', levels(clu), value = T), grep('^Enteroendocrine Cells', levels(clu), value = T), grep('^Tuft Cells', levels(clu), value = T), grep('^Fibro', levels(clu), value = T), grep('^Immune', levels(clu), value = T))
DefaultAssay(clu) <- 'peaks'
annos <- Annotation(clu)
annos <- as.data.frame(annos)
total_peaks <- rownames(clu)
test <- StringToGRanges(total_peaks)
library(BSgenome.Mmusculus.ensembl.mm39)
mouse_views2 <- BSgenomeViews(BSgenome.Mmusculus.ensembl.mm39, test)
#clu$smad <- 0
#clu$tead <- 0
#clu$smad_tead <- 0
#clu$ap1 <- 0
#clu$ap1_smad_tead <- 0
#clu$ap1_smadORtead <- 0
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
top_peaks <- rownames(peaks_mat)
query <- rownames(peaks_mat)
df <- StringToGRanges(query)
mouse_motifs2 <- matchMotifs(pwms = pfm_query, subject = df, out = c('matches'), genome = BSgenome.Mmusculus.ensembl.mm39, bg = colSums(letterFrequency(mouse_views2, letters = c('A', 'C', 'G', 'T'), as.prob = F)) / sum(colSums(letterFrequency(mouse_views2, letters = c('A', 'C', 'G', 'T'), as.prob = F))))
matches_per_peak <- as.matrix(motifMatches(mouse_motifs2))
matches_per_peak <- as.data.frame(matches_per_peak)
matches_per_peak$TEAD <- ifelse(matches_per_peak$TEAD3 == T | matches_per_peak$TEAD2 == T | matches_per_peak$TEAD1 == T | matches_per_peak$TEAD4 == T, T, F)
matches_per_peak$smad <- ifelse(matches_per_peak$SMAD2 == T, T, F)
teads_only <- matches_per_peak[matches_per_peak$TEAD == T & matches_per_peak$smad == F,]
smads_only <- matches_per_peak[matches_per_peak$TEAD == F & matches_per_peak$smad == T,]
teads_smads <- matches_per_peak[matches_per_peak$TEAD == T & matches_per_peak$smad == T,]
tead_peaks <- top_peaks[as.numeric(rownames(teads_only))]
smad_peaks <- top_peaks[as.numeric(rownames(smads_only))]
smad_tead_peaks <- top_peaks[as.numeric(rownames(teads_smads))]
clu$smad <- apply(peaks_mat[smad_peaks, ], 2, mean)
clu$tead <- apply(peaks_mat[tead_peaks, ], 2, mean)
clu$smad_tead <- apply(peaks_mat[smad_tead_peaks,], 2, mean)
library(ggplot2)
meta <- clu@meta.data
meta <- meta[,c('new_collapse_idents', 'smad', 'tead', 'smad_tead')]
compdf2 <- data.frame(clusters = unique(meta$new_collapse_idents), ratio = NA, pvalue = NA, type = 'tead')
for(cluster in compdf2$clusters) {
#compdf2[compdf2$clusters == cluster, 'ratio'] = (mean(meta[WhichCells(clu, idents = 'RevSC1-W13'), 'tead'])) / mean(meta[WhichCells(clu, idents = cluster), 'tead'])
compdf2[compdf2$clusters == cluster, 'pvalue'] = t.test(meta[WhichCells(clu, idents = 'RevSC1-W13'), 'tead'], meta[WhichCells(clu, idents = cluster), 'tead'], alternative = 'two.sided', var.equal = T)$p.value
}
compdf3 <- data.frame(clusters = unique(meta$new_collapse_idents), ratio = NA, pvalue = NA, type = 'smad')
for(cluster in compdf3$clusters) {
#compdf3[compdf3$clusters == cluster, 'ratio'] = (mean(meta[WhichCells(clu, idents = 'RevSC1-W13'), 'smad'])) / mean(meta[WhichCells(clu, idents = cluster), 'smad'])
compdf3[compdf3$clusters == cluster, 'pvalue'] = t.test(meta[WhichCells(clu, idents = 'RevSC1-W13'), 'smad'], meta[WhichCells(clu, idents = cluster), 'smad'], alternative = 'two.sided', var.equal = T)$p.value
}
compdf4 <- data.frame(clusters = unique(meta$new_collapse_idents), ratio = NA, pvalue = NA, type = 'smad_tead')
for(cluster in compdf4$clusters) {
#compdf4[compdf4$clusters == cluster, 'ratio'] = (mean(meta[WhichCells(clu, idents = 'RevSC1-W13'), 'smad_tead'])) / mean(meta[WhichCells(clu, idents = cluster), 'smad_tead'])
compdf4[compdf4$clusters == cluster, 'pvalue'] = t.test(meta[WhichCells(clu, idents = 'RevSC1-W13'), 'smad_tead'], meta[WhichCells(clu, idents = cluster), 'smad_tead'], alternative = 'two.sided',  var.equal = T)$p.value
}
write.csv(compdf2, file = 'tead_pvalues_clu_0_2_dpi_dblt_filtd_apr25_2023.csv')
write.csv(compdf3, file = 'smad_pvalues_clu_0_2_dpi_dblt_filtd_apr25_2023.csv')
write.csv(compdf4, file = 'smad_tead_pvalues_clu_0_2_dpi_dblt_filtd_apr25_2023.csv')
clu <- RenameIdents(clu, 'RevSC1-W13' = 'RevSC')
VlnPlot(clu, features = 'smad', idents = c('RevSC', 'Crypt Cells', 'CVJ', 'Villus Bottom', 'Villus Middle', 'Villus Top', 'Paneth Cells', 'Goblet Cells', 'Enteroendocrine Cells')) + stat_summary(fun = mean, geom = 'point', size = 25, color = 'skyblue', shape = 95, show.legend = F) + scale_y_continuous(expand = c(0,0), limits = c(0, 0.6))
VlnPlot(clu, features = 'smad', idents = c('RevSC', 'Crypt Cells', 'CVJ', 'Villus Bottom', 'Villus Middle', 'Villus Top', 'Paneth Cells', 'Goblet Cells', 'Enteroendocrine Cells')) + stat_summary(fun = mean, geom = 'point', size = 25, color = 'skyblue', shape = 95, show.legend = F)
VlnPlot(clu, features = 'smad', idents = c('RevSC', 'Crypt Cells', 'CVJ', 'Villus Bottom', 'Villus Middle', 'Villus Top', 'Paneth Cells', 'Goblet Cells', 'Enteroendocrine Cells')) + stat_summary(fun = mean, geom = 'point', size = 25, color = 'skyblue', shape = 95, show.legend = F) + scale_y_continuous(expand = c(0,0), limits = c(0, 0.65))
VlnPlot(clu, features = 'tead', idents = c('RevSC', 'Crypt Cells', 'CVJ', 'Villus Bottom', 'Villus Middle','Villus Top', 'Paneth Cells', 'Goblet Cells', 'Enteroendocrine Cells')) + stat_summary(fun = mean, geom = 'point', size = 25, color = 'skyblue', shape = 95, show.legend = F) + scale_y_continuous(expand = c(0,0), limits = c(0, 0.6))
VlnPlot(clu, features = 'smad_tead', idents = c('RevSC', 'Crypt Cells', 'CVJ', 'Villus Bottom','Villus Middle', 'Villus Top', 'Paneth Cells', 'Goblet Cells', 'Enteroendocrine Cells')) + stat_summary(fun = mean, geom = 'point', size = 25, color = 'skyblue', shape = 95, show.legend = F)+ scale_y_continuous(expand = c(0,0), limits = c(0, 0.8))
savehistory("/media/shyam/external/multiome_clu_2/multi5/clu_0_2_dpi_dblt_filtd_smad_tead_violin_analysis.R")
