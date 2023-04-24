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
if(i %in% c('TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'SMAD2', 'FOS::JUN')){
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
teads_only <- matches_per_peak[matches_per_peak$TEAD == T & matches_per_peak$smad == F & matches_per_peak$`FOS::JUN` ==F,]
smads_only <- matches_per_peak[matches_per_peak$TEAD == F & matches_per_peak$smad == T & matches_per_peak$`FOS::JUN` ==F,]
ap1_only <- matches_per_peak[matches_per_peak$TEAD == F & matches_per_peak$smad == F & matches_per_peak$`FOS::JUN` ==T,]
teads_smads <- matches_per_peak[matches_per_peak$TEAD == T & matches_per_peak$smad == T & matches_per_peak$`FOS::JUN` ==F,]
smad_ap1 <- matches_per_peak[matches_per_peak$TEAD == F & matches_per_peak$smad == T & matches_per_peak$`FOS::JUN` ==T,]
tead_ap1 <- matches_per_peak[matches_per_peak$TEAD == T & matches_per_peak$smad == F & matches_per_peak$`FOS::JUN` ==T,]
smad_tead_ap1 <- matches_per_peak[matches_per_peak$TEAD == T & matches_per_peak$smad == T & matches_per_peak$`FOS::JUN` ==T,]
smad_or_tead_ap1 <- matches_per_peak[(matches_per_peak$TEAD == T | matches_per_peak$smad == T) & matches_per_peak$`FOS::JUN` ==T,]
tead_peaks <- top_peaks[as.numeric(rownames(teads_only))]
smad_peaks <- top_peaks[as.numeric(rownames(smads_only))]
ap1_peaks <- top_peaks[as.numeric(rownames(ap1_only))]
smad_tead_peaks <- top_peaks[as.numeric(rownames(teads_smads))]
smad_ap1_peaks <- top_peaks[as.numeric(rownames(smad_ap1))]
tead_ap1_peaks <- top_peaks[as.numeric(rownames(tead_ap1))]
smad_tead_ap1_peaks <- top_peaks[as.numeric(rownames(smad_tead_ap1))]
smad_or_tead_ap1_peaks <- top_peaks[as.numeric(rownames(smad_or_tead_ap1))]
clu$smad <- apply(peaks_mat[smad_peaks, ], 2, mean)
clu$tead <- apply(peaks_mat[tead_peaks, ], 2, mean)
clu$smad_tead <- apply(peaks_mat[smad_tead_peaks,], 2, mean)
clu$ap1 <- apply(peaks_mat[ap1_peaks,], 2, mean)
clu$ap1_smad <- apply(peaks_mat[smad_ap1_peaks,], 2, mean)
clu$ap1_tead <- apply(peaks_mat[tead_ap1_peaks,], 2, mean)
clu$ap1_smad_tead <- apply(peaks_mat[smad_tead_ap1_peaks,], 2, mean)
clu$smad_or_tead_ap1 <- apply(peaks_mat[smad_or_tead_ap1_peaks,], 2, mean)
library(ggplot2)
meta <- clu@meta.data
meta <- meta[,c('new_collapse_idents', 'smad', 'tead', 'smad_tead', 'ap1', 'ap1_smad', 'ap1_tead', 'ap1_smad_tead', 'smad_or_tead_ap1')]
compdf <- data.frame(clusters = unique(meta$new_collapse_idents), mean_diff = NA, type = 'smad_or_tead_ap1')
for(cluster in compdf$clusters) {
compdf[compdf$clusters == cluster, 'mean_diff'] = (mean(meta[WhichCells(clu, idents = 'RevSC1-W13'), 'smad_or_tead_ap1'])) - mean(meta[WhichCells(clu, idents = cluster), 'smad_or_tead_ap1'])
#compdf[compdf$clusters == cluster, 'pvalue'] = t.test(meta[WhichCells(clu, idents = 'RevSC1-W13'), 'smad_or_tead_ap1'], meta[WhichCells(clu, idents = cluster), 'smad_or_tead_ap1'], alternative = 'two.sided',  var.equal = T)$p.value
}
compdf2 <- data.frame(clusters = unique(meta$new_collapse_idents), mean_diff = NA, type = 'smad')
for(cluster in compdf2$clusters) {
compdf2[compdf2$clusters == cluster, 'mean_diff'] = (mean(meta[WhichCells(clu, idents = 'RevSC1-W13'), 'smad'])) - mean(meta[WhichCells(clu, idents = cluster), 'smad'])
#compdf[compdf$clusters == cluster, 'pvalue'] = t.test(meta[WhichCells(clu, idents = 'RevSC1-W13'), 'smad_or_tead_ap1'], meta[WhichCells(clu, idents = cluster), 'smad_or_tead_ap1'], alternative = 'two.sided',  var.equal = T)$p.value
}
compdf3 <- data.frame(clusters = unique(meta$new_collapse_idents), mean_diff = NA, type = 'smad_tead')
for(cluster in compdf3$clusters) {
compdf3[compdf3$clusters == cluster, 'mean_diff'] = (mean(meta[WhichCells(clu, idents = 'RevSC1-W13'), 'smad_tead'])) - mean(meta[WhichCells(clu, idents = cluster), 'smad_tead'])
#compdf[compdf$clusters == cluster, 'pvalue'] = t.test(meta[WhichCells(clu, idents = 'RevSC1-W13'), 'smad_or_tead_ap1'], meta[WhichCells(clu, idents = cluster), 'smad_or_tead_ap1'], alternative = 'two.sided',  var.equal = T)$p.value
}
compdf4 <- data.frame(clusters = unique(meta$new_collapse_idents), mean_diff = NA, type = 'tead')
for(cluster in compdf4$clusters) {
compdf4[compdf4$clusters == cluster, 'mean_diff'] = (mean(meta[WhichCells(clu, idents = 'RevSC1-W13'), 'tead'])) - mean(meta[WhichCells(clu, idents = cluster), 'tead'])
#compdf[compdf$clusters == cluster, 'pvalue'] = t.test(meta[WhichCells(clu, idents = 'RevSC1-W13'), 'smad_or_tead_ap1'], meta[WhichCells(clu, idents = cluster), 'smad_or_tead_ap1'], alternative = 'two.sided',  var.equal = T)$p.value
}
meandf <- data.frame(clusters = unique(meta$new_collapse_idents), mean = NA, var = NA, type = 'smad_or_tead_ap1')
for(cluster in meandf$clusters){
meandf[meandf$clusters == cluster, 'mean'] = mean(meta[WhichCells(clu, idents = cluster), 'smad_or_tead_ap1'])
meandf[meandf$clusters == cluster, 'var'] = var(meta[WhichCells(clu, idents = cluster), 'smad_or_tead_ap1'])
}
meandf2 <- data.frame(clusters = unique(meta$new_collapse_idents), mean = NA, var = NA, type = 'smad')
for(cluster in meandf2$clusters){
meandf2[meandf2$clusters == cluster, 'mean'] = mean(meta[WhichCells(clu, idents = cluster), 'smad'])
meandf2[meandf2$clusters == cluster, 'var'] = var(meta[WhichCells(clu, idents = cluster), 'smad'])
}
meandf3 <- data.frame(clusters = unique(meta$new_collapse_idents), mean = NA, var = NA, type = 'tead')
for(cluster in meandf3$clusters){
meandf3[meandf3$clusters == cluster, 'mean'] = mean(meta[WhichCells(clu, idents = cluster), 'tead'])
meandf3[meandf3$clusters == cluster, 'var'] = var(meta[WhichCells(clu, idents = cluster), 'tead'])
}
meandf4 <- data.frame(clusters = unique(meta$new_collapse_idents), mean = NA, var = NA, type = 'smad_tead')
for(cluster in meandf3$clusters){
meandf4[meandf4$clusters == cluster, 'mean'] = mean(meta[WhichCells(clu, idents = cluster), 'smad_tead'])
meandf4[meandf4$clusters == cluster, 'var'] = var(meta[WhichCells(clu, idents = cluster), 'smad_tead'])
}
library(dplyr)
meandf_main <- bind_rows(meandf, meandf2, meandf3, meandf4)
df <- bind_rows(compdf, compdf2, compdf3, compdf4)
df$type <- factor(df$type, c('smad', 'tead','smad_tead', 'smad_or_tead_ap1'))
df$clusters <- factor(df$clusters, c(unique(grep('RevSC',df$clusters, value = T)), unique(grep('Crypt',df$clusters, value = T)), unique(grep('Paneth',df$clusters, value = T)), unique(grep('CVJ',df$clusters, value = T)), unique(grep('Villus Bottom',df$clusters, value = T)), unique(grep('Villus Middle',df$clusters, value = T)), unique(grep('Villus Top',df$clusters, value = T)), unique(grep('Enteroendocrine Cells',df$clusters, value = T)), unique(grep('Goblet',df$clusters, value = T)), unique(grep('Tuft Cells',df$clusters, value = T)), unique(grep('Fibro',df$clusters, value = T)), unique(grep('Immune',df$clusters, value = T))))
df2 <- df[-c(grep('RevSC', df$clusters), grep('Tuft', df$clusters), grep('Immune', df$clusters), grep('Fibro', df$clusters)),]
df2$var <- NA
for(cluster in unique(df2$clusters)){
for(comp in unique(df2$type)){
df2[df2$clusters == cluster & df2$type == comp, 'var'] = meandf_main[meandf_main$clusters == 'RevSC1-W13' & meandf_main$type == comp, 'var'] + meandf_main[meandf_main$clusters == cluster & meandf_main$type == comp, 'var']
}
}
df3 <- df2
df3$pval <- NA
df3$var <- NA
for(cluster in unique(df3$clusters)){
for(comp in unique(df3$type)){
xstats <- as.data.frame(summary(manova(cbind(smad_tead, tead, smad, smad_or_tead_ap1) ~ new_collapse_idents, data = meta[WhichCells(clu, idents = c('RevSC1-W13', cluster)),]))$stats)
df3[df3$clusters == cluster & df3$type == comp, 'mean_diff'] = df2[df2$clusters == cluster & df2$type == comp, 'mean_diff'] / df2[df2$clusters == cluster & df2$type == 'smad', 'mean_diff']
df3[df3$clusters == cluster & df3$type == comp, 'pval'] = xstats$`Pr(>F)`[1]
df3[df3$clusters == cluster & df3$type == comp, 'var'] = (meandf_main[meandf_main$clusters == 'RevSC1-W13' & meandf_main$type == comp, 'var'] + meandf_main[meandf_main$clusters == cluster & meandf_main$type == comp, 'var']) /  df2[df2$clusters == cluster & df2$type == 'smad', 'mean_diff']
}
}
ggplot(df3, aes(fill = type, y = mean_diff, x = clusters)) + geom_bar(position= 'dodge', stat = 'identity') + geom_errorbar(aes(ymin = mean_diff - var, ymax = mean_diff+var, group = type), width = 0.2, position = position_dodge(0.8))+ scale_y_continuous(expand = c(0,0), limits = c(0, 3))  & theme(panel.background = element_rect(fill = 'white', colour = NULL), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggplot(df3, aes(fill = type, y = mean_diff, x = clusters)) + geom_bar(position= 'dodge', stat = 'identity') + geom_errorbar(aes(ymin = mean_diff - var, ymax = mean_diff+var, group = type), width = 0.2, position = position_dodge(0.8))+ scale_y_continuous(expand = c(0,0), limits = c(0, 3.5))  & theme(panel.background = element_rect(fill = 'white', colour = NULL), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
write.csv(df3, file = 'updated2_df3_smad_normalized_mean_differences_revsc_clusterX_pvalues.csv')
savehistory("/media/shyam/external/multiome_clu_2/multi5/updated2_df3_smad_normalized_mean_differences_revsc_clusterX_pvalues.R")
