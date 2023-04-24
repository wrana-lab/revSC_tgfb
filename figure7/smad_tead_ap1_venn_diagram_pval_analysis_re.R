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
if(i %in% c('TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'SMAD2', 'FOS::JUN')){
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
teads_only <- matches_per_peak[matches_per_peak$TEAD == T & matches_per_peak$smad == F & matches_per_peak$`FOS::JUN` ==F,]
smads_only <- matches_per_peak[matches_per_peak$TEAD == F & matches_per_peak$smad == T & matches_per_peak$`FOS::JUN` ==F,]
ap1_only <- matches_per_peak[matches_per_peak$TEAD == F & matches_per_peak$smad == F & matches_per_peak$`FOS::JUN` ==T,]
teads_smads <- matches_per_peak[matches_per_peak$TEAD == T & matches_per_peak$smad == T & matches_per_peak$`FOS::JUN` ==F,]
smad_ap1 <- matches_per_peak[matches_per_peak$TEAD == F & matches_per_peak$smad == T & matches_per_peak$`FOS::JUN` ==T,]
tead_ap1 <- matches_per_peak[matches_per_peak$TEAD == T & matches_per_peak$smad == F & matches_per_peak$`FOS::JUN` ==T,]
smad_tead_ap1 <- matches_per_peak[matches_per_peak$TEAD == T & matches_per_peak$smad == T & matches_per_peak$`FOS::JUN` ==T,]
smad_or_tead_ap1 <- matches_per_peak[(matches_per_peak$TEAD == T | matches_per_peak$smad == T) & matches_per_peak$`FOS::JUN` ==T,]
tead_peaks <- top_peaks_revsc[as.numeric(rownames(teads_only))]
smad_peaks <- top_peaks_revsc[as.numeric(rownames(smads_only))]
ap1_peaks <- top_peaks_revsc[as.numeric(rownames(ap1_only))]
smad_tead_peaks <- top_peaks_revsc[as.numeric(rownames(teads_smads))]
smad_ap1_peaks <- top_peaks_revsc[as.numeric(rownames(smad_ap1))]
tead_ap1_peaks <- top_peaks_revsc[as.numeric(rownames(tead_ap1))]
smad_tead_ap1_peaks <- top_peaks_revsc[as.numeric(rownames(smad_tead_ap1))]
smad_or_tead_ap1_peaks <- top_peaks_revsc[as.numeric(rownames(smad_or_tead_ap1))]
library(VennDiagram)
xgames <- list(teads = rownames(matches_per_peak[matches_per_peak$TEAD == T,]), smads= rownames(matches_per_peak[matches_per_peak$smad == T,]), ap1 = rownames(matches_per_peak[matches_per_peak$`FOS::JUN` == T,]))
venn.diagram(xgames, filename = 'venn_smad_tead_ap1_redo.tiff')
df_total <- rownames(peaks_mat)
df_total <- StringToGRanges(df_total)
mouse_motifstotal <- matchMotifs(pwms = pfm_query, subject = df_total, out = c('matches'), genome = BSgenome.Mmusculus.ensembl.mm39, bg = colSums(letterFrequency(mouse_views2, letters = c('A', 'C', 'G', 'T'), as.prob = F)) / sum(colSums(letterFrequency(mouse_views2, letters = c('A', 'C', 'G', 'T'), as.prob = F))))
#matches per every peak in global peakset
matches_per_peak_total <- as.matrix(motifMatches(mouse_motifstotal))
matches_per_peak_total <- as.data.frame(matches_per_peak_total)
matches_per_peak_total$TEAD <- ifelse(matches_per_peak_total$TEAD3 == T | matches_per_peak_total$TEAD2 == T | matches_per_peak_total$TEAD1 == T | matches_per_peak_total$TEAD4 == T, T, F)
matches_per_peak_total$smad <- ifelse( matches_per_peak_total$SMAD2 == T, T, F)
smads_only_total <- rownames(matches_per_peak_total[matches_per_peak_total$TEAD == F & matches_per_peak_total$smad == T & matches_per_peak_total$`FOS::JUN` == F ,])
#peaks with smad only
m = length(smads_only_total)
#peaks without smad only
n = nrow(matches_per_peak_total) - m
#diff peaks with smad only
x = c(0:m)
k = length(top_peaks_revsc)
smad_only_probs <- dhyper(x, m, n, k, log = FALSE)
smad_pval <- sum(smad_only_probs[(length(smad_peaks)+1):(m+1)])
smad_pval
smad_pval
#smad_pval is 0.01258392
tead_only_total <- rownames(matches_per_peak_total[matches_per_peak_total$TEAD == T & matches_per_peak_total$smad == F & matches_per_peak_total$`FOS::JUN` == F,])
#peaks with smad only
m = length(tead_only_total)
#peaks without smad only
n = nrow(matches_per_peak_total) - m
#diff peaks with smad only
x = c(0:m)
k = length(top_peaks_revsc)
tead_only_probs <- dhyper(x, m, n, k, log = FALSE)
tead_pval <- sum(tead_only_probs[(length(tead_peaks)+1):(m+1)])
tead_pval
#tead_pval is 3.070918e-06
#ap1 onlt
ap1_only_total <- rownames(matches_per_peak_total[matches_per_peak_total$TEAD == F & matches_per_peak_total$smad == F & matches_per_peak_total$`FOS::JUN` == T,])
#peaks with smad only
m = length(ap1_only_total)
#peaks without smad only
n = nrow(matches_per_peak_total) - m
#diff peaks with smad only
x = c(0:m)
k = length(top_peaks_revsc)
ap1_only_probs <- dhyper(x, m, n, k, log = FALSE)
ap1_pval <- sum(ap1_only_probs[(length(ap1_peaks)+1):(m+1)])
ap1_pval
#ap1 pval is 1.705303e-146
smad_tead_total <- rownames(matches_per_peak_total[matches_per_peak_total$TEAD == T & matches_per_peak_total$smad == T & matches_per_peak_total$`FOS::JUN` == F,])
#peaks with smad only
m = length(smad_tead_total)
#peaks without smad only
n = nrow(matches_per_peak_total) - m
#diff peaks with smad only
x = c(0:m)
k = length(top_peaks_revsc)
smad_tead_probs <- dhyper(x, m, n, k, log = FALSE)
smad_tead_pval <- sum(smad_tead_probs[(length(smad_tead_peaks)+1):(m+1)])
smad_tead_pval
#smad_tead pval is 0.0003617299
#smad_ap1
smad_ap1_total <- rownames(matches_per_peak_total[matches_per_peak_total$TEAD == F & matches_per_peak_total$smad == T & matches_per_peak_total$`FOS::JUN` == T,])
#peaks with smad only
m = length(smad_ap1_total)
#peaks without smad only
n = nrow(matches_per_peak_total) - m
#diff peaks with smad only
x = c(0:m)
k = length(top_peaks_revsc)
smad_ap1_probs <- dhyper(x, m, n, k, log = FALSE)
smad_ap1_pval <- sum(smad_ap1_probs[(length(smad_ap1_peaks)+1):(m+1)])
smad_ap1_pval
smad_ap1_pval
#smad_ap1_pval is 1.715599e-41
#tead_ap1
tead_ap1_total <- rownames(matches_per_peak_total[matches_per_peak_total$TEAD == T & matches_per_peak_total$smad == F & matches_per_peak_total$`FOS::JUN` == T,])
#peaks with smad only
m = length(tead_ap1_total)
#peaks without smad only
n = nrow(matches_per_peak_total) - m
#diff peaks with smad only
x = c(0:m)
k = length(top_peaks_revsc)
tead_ap1_probs <- dhyper(x, m, n, k, log = FALSE)
tead_ap1_pval <- sum(tead_ap1_probs[(length(tead_ap1_peaks)+1):(m+1)])
tead_ap1_pval
#tead_ap1_pval is 9.263116e-76
#tead_ap1_smad
smad_tead_ap1_total <- rownames(matches_per_peak_total[matches_per_peak_total$TEAD == T & matches_per_peak_total$smad == T & matches_per_peak_total$`FOS::JUN` == T,])
#peaks with smad only
m = length(smad_tead_ap1_total)
#peaks without smad only
n = nrow(matches_per_peak_total) - m
#diff peaks with smad only
x = c(0:m)
k = length(top_peaks_revsc)
smad_tead_ap1_probs <- dhyper(x, m, n, k, log = FALSE)
smad_tead_ap1_pval <- sum(smad_tead_ap1_probs[(length(smad_tead_ap1_peaks)+1):(m+1)])
smad_tead_ap1_pval
#smad_tead_ap1_pval is 7.549338e-12
savehistory("/media/shyam/external/multiome_clu_2/multi5/smad_tead_ap1_venn_diagram_pval_analysis_re.R")
