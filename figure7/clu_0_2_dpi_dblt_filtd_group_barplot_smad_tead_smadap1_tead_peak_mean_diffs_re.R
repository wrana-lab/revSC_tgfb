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
metaLM_smad <- lm(smad ~ new_collapse_idents, data = meta)
smad_anova <- anova(metaLM_smad)
smad_tukey <- TukeyHSD(aov(metaLM_smad))
metaLM_tead<- lm(tead ~ new_collapse_idents, data = meta)
tead_anova <- anova(metaLM_tead)
tead_tukey <- TukeyHSD(aov(metaLM_tead))
metaLM_smad_tead<- lm(smad_tead ~ new_collapse_idents, data = meta)
smad_tead_anova <- anova(metaLM_smad_tead)
smad_tead_tukey <- TukeyHSD(aov(metaLM_smad_tead))
metaLM_ap1<- lm(ap1 ~ new_collapse_idents, data = meta)
ap1_anova <- anova(metaLM_ap1)
ap1_tukey <- TukeyHSD(aov(metaLM_ap1))
metaLM_ap1_smad <- lm(ap1_smad ~ new_collapse_idents, data = meta)
ap1_smad_anova <- anova(metaLM_ap1_smad)
ap1_smad_tukey <- TukeyHSD(aov(metaLM_ap1_smad))
metaLM_ap1_tead <- lm(ap1_tead ~ new_collapse_idents, data = meta)
ap1_tead_anova <- anova(metaLM_ap1_tead)
ap1_tead_tukey <- TukeyHSD(aov(metaLM_ap1_tead))
metaLM_ap1_smad_tead <- lm(ap1_smad_tead ~ new_collapse_idents, data = meta)
ap1_smad_tead_anova <- anova(metaLM_ap1_smad_tead)
ap1_smad_tead_tukey <- TukeyHSD(aov(metaLM_ap1_smad_tead))
metaLM_smad_or_tead_ap1 <- lm(smad_or_tead_ap1 ~ new_collapse_idents, data = meta)
smad_or_tead_ap1_anova <- anova(metaLM_smad_or_tead_ap1)
smad_or_tead_ap1_tukey <- TukeyHSD(aov(metaLM_smad_or_tead_ap1))
smad_tukeyi <- smad_tukey$new_collapse_idents[unique(c(grep('RevSC',rownames(smad_tukey$new_collapse_idents)), grep('Crypt',rownames(smad_tukey$new_collapse_idents)))),]
tead_tukeyi <- tead_tukey$new_collapse_idents[unique(c(grep('RevSC',rownames(tead_tukey$new_collapse_idents)), grep('Crypt',rownames(tead_tukey$new_collapse_idents)))),]
smad_tead_tukeyi <- smad_tead_tukey$new_collapse_idents[unique(c(grep('RevSC',rownames(smad_tead_tukey$new_collapse_idents)), grep('Crypt',rownames(smad_tead_tukey$new_collapse_idents)))),]
ap1_tukeyi <- ap1_tukey$new_collapse_idents[unique(c(grep('RevSC',rownames(ap1_tukey$new_collapse_idents)), grep('Crypt',rownames(ap1_tukey$new_collapse_idents)))),]
ap1_smad_tukeyi <- ap1_smad_tukey$new_collapse_idents[unique(c(grep('RevSC',rownames(ap1_smad_tukey$new_collapse_idents)), grep('Crypt',rownames(ap1_smad_tukey$new_collapse_idents)))),]
ap1_tead_tukeyi <- ap1_tead_tukey$new_collapse_idents[unique(c(grep('RevSC',rownames(ap1_tead_tukey$new_collapse_idents)), grep('Crypt',rownames(ap1_tead_tukey$new_collapse_idents)))),]
ap1_smad_tead_tukeyi <- ap1_smad_tead_tukey$new_collapse_idents[unique(c(grep('RevSC',rownames(ap1_smad_tead_tukey$new_collapse_idents)), grep('Crypt',rownames(ap1_smad_tead_tukey$new_collapse_idents)))),]
smad_or_tead_ap1_tukeyi <- smad_or_tead_ap1_tukey$new_collapse_idents[unique(c(grep('RevSC',rownames(smad_or_tead_ap1_tukey$new_collapse_idents)), grep('Crypt',rownames(smad_or_tead_ap1_tukey$new_collapse_idents)))),]
smad_or_tead_ap1_tukeyi
smad_tukeyi <- as.data.frame(smad_tukey$new_collapse_idents[unique(c(grep('RevSC',rownames(smad_tukey$new_collapse_idents)))),])
tead_tukeyi <- as.data.frame(tead_tukey$new_collapse_idents[unique(c(grep('RevSC',rownames(tead_tukey$new_collapse_idents)))),])
smad_tead_tukeyi <- as.data.frame(smad_tead_tukey$new_collapse_idents[unique(c(grep('RevSC',rownames(smad_tead_tukey$new_collapse_idents)))),])
ap1_tukeyi <- as.data.frame(ap1_tukey$new_collapse_idents[unique(c(grep('RevSC',rownames(ap1_tukey$new_collapse_idents)))),])
ap1_smad_tukeyi <- as.data.frame(ap1_smad_tukey$new_collapse_idents[unique(c(grep('RevSC',rownames(ap1_smad_tukey$new_collapse_idents)))),])
ap1_tead_tukeyi <- as.data.frame(ap1_tead_tukey$new_collapse_idents[unique(c(grep('RevSC',rownames(ap1_tead_tukey$new_collapse_idents)))),])
ap1_smad_tead_tukeyi <- as.data.frame(ap1_smad_tead_tukey$new_collapse_idents[unique(c(grep('RevSC',rownames(ap1_smad_tead_tukey$new_collapse_idents)))),])
smad_or_tead_ap1_tukeyi <- as.data.frame(smad_or_tead_ap1_tukey$new_collapse_idents[unique(c(grep('RevSC',rownames(smad_or_tead_ap1_tukey$new_collapse_idents)))),])
smad_tukeyi$group <- 'smad'
tead_tukeyi$group <- 'tead'
smad_tead_tukeyi$group <- 'smad_tead'
ap1_tukeyi$group <- 'ap1'
ap1_smad_tukeyi$group <- 'ap1_smad'
ap1_tead_tukeyi$group <- 'ap1_tead'
ap1_smad_tead_tukeyi$group <- 'ap1_smad_tead'
smad_or_tead_ap1_tukeyi$group <- 'smad_or_tead_ap1'
smad_tukeyi$clusters <- rownames(smad_tukeyi)
tead_tukeyi$clusters <- rownames(tead_tukeyi)
smad_tead_tukeyi$clusters <- rownames(smad_tead_tukeyi)
ap1_tukeyi$clusters <- rownames(ap1_tukeyi)
ap1_smad_tukeyi$clusters <- rownames(ap1_smad_tukeyi)
ap1_tead_tukeyi$clusters <- rownames(ap1_tead_tukeyi)
ap1_smad_tead_tukeyi$clusters <- rownames(ap1_smad_tead_tukeyi)
smad_or_tead_ap1_tukeyi$clusters <- rownames(smad_or_tead_ap1_tukeyi)
rownames(smad_tukeyi) <- 1:nrow(smad_tukeyi)
rownames(tead_tukeyi) <- 1:nrow(tead_tukeyi)
rownames(smad_tead_tukeyi) <- 1:nrow(smad_tead_tukeyi)
rownames(ap1_tukeyi) <- 1:nrow(ap1_tukeyi)
rownames(ap1_smad_tukeyi) <- 1:nrow(ap1_smad_tukeyi)
rownames(ap1_tead_tukeyi) <- 1:nrow(ap1_tead_tukeyi)
rownames(ap1_smad_tead_tukeyi) <- 1:nrow(ap1_smad_tead_tukeyi)
rownames(smad_or_tead_ap1_tukeyi) <- 1:nrow(smad_or_tead_ap1_tukeyi)
library(dplyr)
df <- bind_rows(smad_tukeyi, tead_tukeyi, smad_or_tead_ap1_tukeyi)
df$lwr <- NULL
df$upr <- NULL
df$diff <- abs(df$diff)
df$group <- factor(df$group, c('smad', 'tead', 'smad_or_tead_ap1'))
df$clusters <- factor(df$clusters, c(unique(grep('Crypt',df$clusters, value = T)), unique(grep('Paneth',df$clusters, value = T)), unique(grep('CVJ',df$clusters, value = T)), unique(grep('Villus Bottom',df$clusters, value = T)), unique(grep('Villus Middle',df$clusters, value = T)), unique(grep('Villus Top',df$clusters, value = T)), unique(grep('Enteroendocrine Cells',df$clusters, value = T)), unique(grep('Goblet',df$clusters, value = T)), unique(grep('Tuft Cells',df$clusters, value = T)), unique(grep('Fibro',df$clusters, value = T)), unique(grep('Immune',df$clusters, value = T))))
ggplot(df, aes(fill = group, y = diff, x = clusters)) + geom_bar(position= 'dodge', stat = 'identity') + scale_y_continuous(expand = c(0,0), limits = c(0, 0.082)) & theme(panel.background = element_rect(fill = 'white', colour = NULL), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
df2 <- df[-c(grep('Tuft', df$clusters), grep('Immune', df$clusters), grep('Fibro', df$clusters)),]
ggplot(df2, aes(fill = group, y = diff, x = clusters)) + geom_bar(position= 'dodge', stat = 'identity') + scale_y_continuous(expand = c(0,0), limits = c(0, 0.078)) & theme(panel.background = element_rect(fill = 'white', colour = NULL), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
levels(clu)
clu <- RenameIdents(clu, 'RevSC1-W13' = 'RevSC')
VlnPlot(clu, features = 'smad', idents = c('RevSC', 'Crypt Cells', 'CVJ', 'Villus Bottom', 'Villus Middle', 'Villus Top', 'Paneth Cells', 'Goblet Cells', 'Enteroendocrine Cells')) + stat_summary(fun = mean, geom = 'point', size = 25, color = 'skyblue', shape = 95, show.legend = F) + scale_y_continuous(expand = c(0,0), limits = c(0, 0.6))
VlnPlot(clu, features = 'tead', idents = c('RevSC', 'Crypt Cells', 'CVJ', 'Villus Bottom', 'Villus Middle','Villus Top', 'Paneth Cells', 'Goblet Cells', 'Enteroendocrine Cells')) + stat_summary(fun = mean, geom = 'point', size = 25, color = 'skyblue', shape = 95, show.legend = F) + scale_y_continuous(expand = c(0,0), limits = c(0, 0.6))
VlnPlot(clu, features = 'smad_tead', idents = c('RevSC', 'Crypt Cells', 'CVJ', 'Villus Bottom','Villus Middle', 'Villus Top', 'Paneth Cells', 'Goblet Cells', 'Enteroendocrine Cells')) + stat_summary(fun = mean, geom = 'point', size = 25, color = 'skyblue', shape = 95, show.legend = F)+ scale_y_continuous(expand = c(0,0), limits = c(0, 0.8))
VlnPlot(clu, features = 'smad_or_tead_ap1', idents = c('RevSC', 'Crypt Cells', 'CVJ', 'Villus Bottom','Villus Middle', 'Villus Top', 'Paneth Cells', 'Goblet Cells', 'Enteroendocrine Cells')) + stat_summary(fun = mean, geom = 'point', size = 25, color = 'skyblue', shape = 95, show.legend = F)+ scale_y_continuous(expand = c(0,0), limits = c(0, 0.8))
savehistory("/media/shyam/external/multiome_clu_2/multi5/clu_0_2_dpi_dblt_filtd_group_barplot_smad_tead_smadap1_tead_peak_mean_diffs_re.R")
