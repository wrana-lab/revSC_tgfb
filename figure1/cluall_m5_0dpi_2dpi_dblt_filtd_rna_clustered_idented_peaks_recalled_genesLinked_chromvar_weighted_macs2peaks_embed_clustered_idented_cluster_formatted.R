setwd("/media/shyam/external/multiome_clu_2/multi5")
library(Seurat)
library(Signac)
load("/media/shyam/external/multiome_clu_2/multi5/cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar_weighted_macs2peaks_embed_clustered_idented.RData")
Idents(clu) <- clu$peaks_clusters_idents
df_peaks <- data.frame(peaksUMAP1 = Embeddings(clu, 'umap.peaks')[,1], peaksUMAP2 = Embeddings(clu, 'umap.peaks')[,2], ClusterNumber = clu$peaks_clusters, ClusterName = clu$peaks_clusters_idents)
df_label_peaks <- df_peaks %>% group_by(ClusterNumber) %>% summarize(peaksUMAP1 = mean(peaksUMAP1), peaksUMAP2 = mean(peaksUMAP2))
library(tidyverse)
df_label_peaks <- df_peaks %>% group_by(ClusterNumber) %>% summarize(peaksUMAP1 = mean(peaksUMAP1), peaksUMAP2 = mean(peaksUMAP2))
ggplot(df_peaks, aes(x = peaksUMAP1, y = peaksUMAP2)) + geom_point(shape = 16, aes(color = ClusterName)) + geom_text(data = df_label_peaks, aes(label = ClusterNumber)) + theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 1, color = 'black'), axis.text = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), text = element_text(size = 12), legend.key = element_rect(fill = 'white')) + guides(color = guide_legend(override.aes = list(size = 3) ) ) + labs(color = 'ATAC Clusters')
ggplot(df_peaks, aes(x = peaksUMAP1, y = peaksUMAP2)) + geom_point(shape = 12, aes(color = ClusterName)) + geom_text(data = df_label_peaks, aes(label = ClusterNumber)) + theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 1, color = 'black'), axis.text = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), text = element_text(size = 12), legend.key = element_rect(fill = 'white')) + guides(color = guide_legend(override.aes = list(size = 3) ) ) + labs(color = 'ATAC Clusters')
?geom_point
ggplot(df_peaks, aes(x = peaksUMAP1, y = peaksUMAP2)) + geom_point(shape = 16, aes(color = ClusterName), size = 12) + geom_text(data = df_label_peaks, aes(label = ClusterNumber)) + theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 1, color = 'black'), axis.text = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), text = element_text(size = 12), legend.key = element_rect(fill = 'white')) + guides(color = guide_legend(override.aes = list(size = 3) ) ) + labs(color = 'ATAC Clusters')
ggplot(df_peaks, aes(x = peaksUMAP1, y = peaksUMAP2)) + geom_point(shape = 16, aes(color = ClusterName), size = 5) + geom_text(data = df_label_peaks, aes(label = ClusterNumber)) + theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 1, color = 'black'), axis.text = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), text = element_text(size = 12), legend.key = element_rect(fill = 'white')) + guides(color = guide_legend(override.aes = list(size = 3) ) ) + labs(color = 'ATAC Clusters')
ggplot(df_peaks, aes(x = peaksUMAP1, y = peaksUMAP2)) + geom_point(shape = 16, aes(color = ClusterName), size = 1) + geom_text(data = df_label_peaks, aes(label = ClusterNumber)) + theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 1, color = 'black'), axis.text = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), text = element_text(size = 12), legend.key = element_rect(fill = 'white')) + guides(color = guide_legend(override.aes = list(size = 3) ) ) + labs(color = 'ATAC Clusters')
ggplot(df_peaks, aes(x = peaksUMAP1, y = peaksUMAP2)) + geom_point(shape = 16, aes(color = ClusterName), size = 0.5) + geom_text(data = df_label_peaks, aes(label = ClusterNumber)) + theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 1, color = 'black'), axis.text = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), text = element_text(size = 12), legend.key = element_rect(fill = 'white')) + guides(color = guide_legend(override.aes = list(size = 3) ) ) + labs(color = 'ATAC Clusters')
ggplot(df_peaks, aes(x = peaksUMAP1, y = peaksUMAP2)) + geom_point(shape = 16, aes(color = ClusterName), size = 1) + geom_text(data = df_label_peaks, aes(label = ClusterNumber)) + theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 1, color = 'black'), axis.text = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), text = element_text(size = 12), legend.key = element_rect(fill = 'white')) + guides(color = guide_legend(override.aes = list(size = 3) ) ) + labs(color = 'ATAC Clusters')
Idents(clu) <- clu$wnn_clusters_idents
for(i in c('EC', 'GC', 'EE', 'PC', 'TC','CC' ,'Immune', 'Fibro')){
clu <- SetIdent(clu, cells = WhichCells(clu, expression = wnn_clusters_idents %in% grep(paste('^', i, '.*', sep = ""), unique(clu$wnn_clusters_idents), value = T)), value = i)
}
Idents(clu) <- clu$wnn_clusters_idents
for(i in c('EC', 'GC', 'EE', 'PC', 'TC','ICC' ,'Immune', 'Fibro')){
clu <- SetIdent(clu, cells = WhichCells(clu, expression = wnn_clusters_idents %in% grep(paste('^', i, '.*', sep = ""), unique(clu$wnn_clusters_idents), value = T)), value = i)
}
DimPlot(clu, reduction = 'umap.rna', label = T)
DimPlot(clu, reduction = 'wnn.umap', label = T)
clu$new_wnn_clusters <- paste(Idents(clu), '-','W', clu$wnn_clusters, sep = "")
df <- data.frame(wnnUMAP1 = Embeddings(clu, 'wnn.umap')[,1], wnnUMAP2 = Embeddings(clu, 'wnn.umap')[,2], ClusterNumber = clu$wnn_clusters, ClusterName = clu$new_wnn_clusters)
df_label <- df %>% group_by(ClusterNumber) %>% summarize(wnnUMAP1 = mean(wnnUMAP1), wnnUMAP2 = mean(wnnUMAP2))
df$ClusterName <- factor(df$ClusterName, levels=levels(clu))
ggplot(df, aes(x = wnnUMAP1, y = wnnUMAP2)) + geom_point(shape = 16, aes(color = ClusterName)) + geom_text(data = df_label, aes(label = ClusterNumber))
ggplot(df, aes(x = wnnUMAP1, y = wnnUMAP2)) + geom_point(shape = 16, aes(color = ClusterName)) + geom_text(data = df_label, aes(label = ClusterNumber)) + theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 1, color = 'black'), axis.text = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), text = element_text(size = 12), legend.key = element_rect(fill = 'white')) + guides(color = guide_legend(override.aes = list(size = 3) ) ) + labs(color = 'WNN Clusters')
View(df)
df <- data.frame(wnnUMAP1 = Embeddings(clu, 'wnn.umap')[,1], wnnUMAP2 = Embeddings(clu, 'wnn.umap')[,2], ClusterNumber = clu$wnn_clusters, ClusterName = clu$new_wnn_clusters)
df_label <- df %>% group_by(ClusterNumber) %>% summarize(wnnUMAP1 = mean(wnnUMAP1), wnnUMAP2 = mean(wnnUMAP2))
ggplot(df, aes(x = wnnUMAP1, y = wnnUMAP2)) + geom_point(shape = 16, aes(color = ClusterName)) + geom_text(data = df_label, aes(label = ClusterNumber)) + theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 1, color = 'black'), axis.text = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), text = element_text(size = 12), legend.key = element_rect(fill = 'white')) + guides(color = guide_legend(override.aes = list(size = 3) ) ) + labs(color = 'WNN Clusters')
ggplot(df, aes(x = wnnUMAP1, y = wnnUMAP2)) + geom_point(shape = 16, aes(color = ClusterName), size = 1) + geom_text(data = df_label, aes(label = ClusterNumber)) + theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 1, color = 'black'), axis.text = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), text = element_text(size = 12), legend.key = element_rect(fill = 'white')) + guides(color = guide_legend(override.aes = list(size = 3) ) ) + labs(color = 'WNN Clusters')
DimPlot(clu, reduction = 'wnn.umap', split.by = 'time') & scale_fill_manual(values = c('red', 'blue'))
DimPlot(clu, reduction = 'wnn.umap', group.by = 'time') & scale_fill_manual(values = c('red', 'blue'))
df_time <- data.frame(wnnUMAP1 = Embeddings(clu, 'wnn.umap')[,1], wnnUMAP2 = Embeddings(clu, 'wnn.umap')[,2], time = clu$time)
ggplot(df_time, aes(x = wnnUMAP1, y = wnnUMAP2)) + geom_point(shape = 16, aes(color = time)) + theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 1, color = 'black'), axis.text = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), text = element_text(size = 12), legend.key = element_rect(fill = 'white')) + guides(color = guide_legend(override.aes = list(size = 3) ) ) + labs(color = 'WNN Clusters')
ggplot(df_time, aes(x = wnnUMAP1, y = wnnUMAP2)) + geom_point(shape = 16, aes(color = time)) + scale_color_manual(values = c('red', 'blue'))+ theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 1, color = 'black'), axis.text = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), text = element_text(size = 12), legend.key = element_rect(fill = 'white')) + guides(color = guide_legend(override.aes = list(size = 3) ) ) + labs(color = 'WNN Clusters')
ggplot(df_time, aes(x = wnnUMAP1, y = wnnUMAP2)) + geom_point(shape = 16, aes(color = time)) + scale_color_manual(values = c('#2166AC', '#B2182B'))+ theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 1, color = 'black'), axis.text = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), text = element_text(size = 12), legend.key = element_rect(fill = 'white')) + guides(color = guide_legend(override.aes = list(size = 3) ) ) + labs(color = 'WNN Clusters')
ggplot(df_time, aes(x = wnnUMAP1, y = wnnUMAP2)) + geom_point(shape = 16, aes(color = time), size = 1) + scale_color_manual(values = c('#2166AC', '#B2182B'))+ theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 1, color = 'black'), axis.text = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), text = element_text(size = 12), legend.key = element_rect(fill = 'white')) + guides(color = guide_legend(override.aes = list(size = 3) ) ) + labs(color = 'WNN Clusters')
ggplot(df_time, aes(x = wnnUMAP1, y = wnnUMAP2)) + geom_point(shape = 16, aes(color = time), size = 1) + scale_color_manual(values = c('#2166AC', '#B2182B'))+ theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 1, color = 'black'), axis.text = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), text = element_text(size = 12), legend.key = element_rect(fill = 'white')) + guides(color = guide_legend(override.aes = list(size = 3) ) ) + labs(color = 'Time')
revsc_markers <- c("Malat1",
"Ly6d",
"Clu",
"Cldn4",
"Sprr1a",
"Areg",
"Gsta1",
"Tm4sf4",
"Lypd8",
"Neat1",
"Ubd",
"Ctsd",
"Reg3g",
"Cdkn1a",
"Guca2a",
"Emp1",
"Muc3",
"Ms4a10",
"Gm42418",
"Ahnak",
"Lamc2",
"Serpinb1a",
"Cdh1",
"Krt19",
"Cdhr5",
"Ccnd2",
"Cldn3",
"Cxadr",
"Lgals3",
"Gsta4",
"F3",
"Cxcl16",
"2210407C18Rik",
"F11r",
"Xist",
"Reg3b",
"Prap1",
"Guca2b",
"Itm2b",
"Nupr1",
"Anxa1",
"Basp1",
"Krt8",
"Tob1",
"Itgb4",
"Cdh17",
"Anxa2",
"Plec",
"Ms4a8a",
"Gsn")
DefaultAssay(clu) <- 'RNA'
top10_revsc <- revsc_markers[revsc_markers %in% rownames(clu)][1:10]
clu$revsc_sig <- apply(clu@assays$RNA@data[top10_revsc,], MARGIN = 2, median)
DefaultAssay(clu) <- 'chromvar'
motiffss <- as.data.frame(t(as.data.frame(clu@assays$peaks@motifs@motif.names)))
rownames(clu@assays$chromvar@data) <- motiffss$V1
Idents(clu) <- clu$new_wnn_clusters
DimPlot(clu, reduction = 'wnn.umap', label = T)
revsc_chrom <- FindMarkers(clu, ident.1 = 'RevSC1-W13', mean.fxn = rowMeans, fc.name = 'avg_diff', only.pos = T)
revsc_chrom$log_pval <-  -log10(revsc_chrom$p_val)
revsc_chrom$ord <- order(revsc_chrom$log_pval)
revsc_chrom$motifs <- rownames(revsc_chrom)
spec_revsc <- revsc_chrom[which(revsc_chrom$motifs %in% c('FOS', 'FOS::JUN', 'JUN(var.2)','SMAD3', 'Smad4', 'Smad2::Smad3', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4')),]
prevsc + geom_point(data = spec_revsc, aes(x = ord, y = log_pval), color = 'red', show.legend = F) + geom_label_repel(data = spec_revsc, aes(label = motifs), box.padding = 0.75, force = 1, segment.color = 'black') + theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 0.5), axis.text = element_text(size = 12)) + xlab('rank order') + ylab('-log10_pval')
prevsc <- ggplot(revsc_chrom, aes(x = ord, y = log_pval,label = motifs)) + geom_point(color = 'darkgray')
prevsc + geom_point(data = spec_revsc, aes(x = ord, y = log_pval), color = 'red', show.legend = F) + geom_label_repel(data = spec_revsc, aes(label = motifs), box.padding = 0.75, force = 1, segment.color = 'black') + theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 0.5), axis.text = element_text(size = 12)) + xlab('rank order') + ylab('-log10_pval')
library(ggrepel)
prevsc + geom_point(data = spec_revsc, aes(x = ord, y = log_pval), color = 'red', show.legend = F) + geom_label_repel(data = spec_revsc, aes(label = motifs), box.padding = 0.75, force = 1, segment.color = 'black') + theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 0.5), axis.text = element_text(size = 12)) + xlab('rank order') + ylab('-log10_pval')
prevsc + geom_point(data = spec_revsc, aes(x = ord, y = log_pval), color = '#B2182B', show.legend = F) + geom_label_repel(data = spec_revsc, aes(label = motifs), box.padding = 0.75, force = 1, segment.color = 'black') + theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(size = 0.5), axis.text = element_text(size = 12)) + xlab('rank order') + ylab('-log10_pval')
DotPlot(clu, features = c('FOS', 'FOS::JUN', 'JUN(var.2)', 'SMAD3', 'Smad4', 'Smad2::Smad3', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4'))
cvj <- c("Car4",
"Ccl25",
"Mgst1",
"Arg2",
"Tstd1",
"Prpsap1",
"Atp5g1",
"Mgst3",
"Ndufb5",
"Cox7c",
"Cox5a",
"Cox7b",
"Cox8a",
"Uqcr11",
"Cox5b",
"Chchd10",
"Atp5c1",
"Atpif1",
"Atp5l",
"Mtch2",
"Ndufa4",
"Ndufb9",
"Mpc2",
"Atp5a1",
"Gna11",
"Prdx1",
"Dbi",
"Sis",
"Maob",
"Slc7a8",
"Gstm3",
"Fam132a",
"Reg3a",
"Arg2",
"Aldh9a1",
"Mgst1",
"Glud1",
"Tstd1",
"Car4",
"Ccl25",
"Prpsap1",
"Cyb5a",
"Prdx1",
"Mgst2",
"Mgst3",
"Gstp1",
"Gna11",
"Sis",
"Reg1",
"Sis",
"Gsta1",
"Gstm3",
"Arg2",
"Ces1f",
"Bche",
"Acaa2",
"2210407C18Rik",
"Fth1",
"Chpt1",
"Il18",
"Chchd10",
"Cyb5a",
"Plac8",
"Mgst3",
"Ccl25",
"Atp5g1",
"Reg3b",
"Reg3g",
"Ftl1",
"Ost4",
"Gng5",
"Plac8",
"Reg3g"
)
villus_bottom <- c("Apol10a",
"Aoc1",
"Rdh7",
"Cyp2b10",
"Apol10a",
"Cyp3a25",
"Slc6a20a",
"Aldh1a1",
"Maoa",
"Aoc1",
"Lct",
"Naaladl1",
"Sult1d1",
"Cox7a1",
"Cat",
"Ano6",
"Ugt2b34",
"Abcd3",
"Dnase1",
"Cndp2",
"Adipor2",
"Adh1",
"Hadh",
"Mogat2",
"Slc43a2",
"Cyp4v3",
"Fabp1",
"Rdh7",
"Cyp2b10",
"Apol10a",
"Cyp3a25",
"Slc6a20a",
"Aldh1a1",
"Maoa",
"Aoc1",
"Lct",
"Naaladl1",
"Sult1d1",
"Cox7a1",
"Cat",
"Ano6",
"Ugt2b34",
"Abcd3",
"Dnase1",
"Cndp2",
"Adipor2",
"Adh1",
"Hadh",
"Mogat2",
"Slc43a2",
"Cyp4v3",
"Fabp1",
"Rdh7",
"Cyp2b10",
"Spink1",
"Tm4sf20",
"Cdh17",
"Ace",
"Prdx5",
"Fbln1",
"Gda",
"Acsl5",
"Dpp4",
"Mgam",
"Asah2",
"Mttp",
"Ces2e",
"Enpep",
"Cyp4v3",
"Fabp2",
"Rbp2",
"Cyp4f14",
"Anpep",
"Leap2",
"Cycs",
"Serf2",
"Rbp2",
"Fabp2",
"Prdx5",
"Abhd11os",
"Smdt1",
"Scp2",
"Gm10116",
"Spink1",
"Calml4",
"Pdcd6",
"Fabp1",
"Gsto1",
"Cda",
"Cript",
"Gda",
"Dhrs1",
"Tpi1",
"Anpep",
"Pls1",
"Slc5a1",
"Mgam",
"Gda",
"Ace",
"Tcn2",
"Leap2",
"Naprt",
"Slc35c2",
"D5Ertd579e",
"Mical1",
"Naaladl1",
"Slc51b",
"Fabp2",
"Enpep",
"Dpp4",
"Anpep",
"Slc51a",
"Abcg2",
"Tm4sf20"
)
villus_mid <- c("St3gal4",
"Ces2a",
"St3gal4",
"Cyp3a11",
"Gapdh",
"Ephx2",
"Clca4b",
"Cyp3a13",
"Ace2",
"Apob",
"Mep1b",
"Prap1",
"Slc15a1",
"Ifit1",
"Guca2b",
"Treh",
"B2m",
"Cst6",
"Muc13",
"Rfk",
"Ifi27l2b",
"H2-Q2",
"Slc6a19",
"2010106E10Rik",
"Sprr2a3",
"Crip1",
"Edf1",
"Rfk",
"Hcfc1r1",
"Max",
"Sprr2a3",
"Actb",
"Ggt1",
"H2-K1",
"Dgat1",
"Prap1",
"Cdhr5",
"Cdhr2",
"H2-Q2",
"Slc6a19",
"Guca2b",
"Gsdmd",
"Slc9a3r1",
"Eps8l2",
"Slc26a6",
"Xdh",
"Cgref1",
"Oit1",
"Adap1",
"Clec2d",
"Ceacam1",
"Clca4b",
"Slc9a3r1",
"Cdhr5",
"App",
"Cdhr2",
"Slc6a19",
"Sepp1",
"Upp1",
"Dhcr24",
"Sprr2a3",
"Slc27a4",
"Apob",
"Slc6a8",
"Muc13",
"Sepp1",
"H2-K1",
"Cdhr2",
"Slc6a19",
"Ceacam1",
"Slc15a1"
)
villus_top <- c("Adh6a",
"Muc3",
"Apoc3",
"Anxa2",
"Krt20",
"Apoa4",
"Apoa1",
"Cystm1",
"2200002D01Rik",
"S100a10",
"Tmsb10",
"Prr13",
"Dstn",
"Gabarap",
"Rhoc",
"Krt20",
"Hist1h2bc",
"Clic1",
"Rac1",
"2010003K11Rik",
"Ada",
"Pxdc1",
"Irf7",
"Myo15b",
"Malat1",
"Apoa4",
"Apoa1",
"2010003K11Rik",
"Map2k2",
"Galnt6",
"Gm42418",
"Pcyt2",
"Abhd2",
"Rnf186",
"Atp10b",
"Pmp22",
"2010109I03Rik",
"Slc34a2",
"Clca4a",
"2200002D01Rik",
"Mxd1",
"Krt20",
"S100a10",
"Prr13",
"Apoa1",
"Ezr",
"Muc3",
"Serpinb1a",
"Npc1l1",
"Fam3b",
"Tm4sf4",
"Slc28a2",
"Clca4a",
"2010109I03Rik",
"Ada",
"2010107G12Rik",
"Apoa4",
"Myo15b",
"Apoc3",
"Apoa1",
"Muc3",
"Mxd1",
"Lmo7",
"Krt20",
"Ceacam20",
"Creb3l3",
"Malat1",
"Noct",
"Plec",
"Sfn",
"Atp10b",
"Abhd2",
"Enpp3",
"Serpinb1a",
"4930539E08Rik",
"Slc34a2",
"Pmp22",
"Klf4",
"Pxdc1",
"Slc25a22",
"Ept1",
"Dnpep",
"Ormdl3",
"Nt5e",
"Lpin2",
"Fosl2",
"Ifrd1"
)
cvj <- unique(cvj)
villus_bottom <- unique(villus_bottom)
villus_mid <- unique(villus_mid)
villus_top <- unique(villus_top)
DefaultAssay(clu) <- 'RNA'
clu$cvj <- apply(clu@assays$RNA@data[cvj[cvj %in% rownames(clu@assays$RNA@data)],], 2, median)
clu$villus_bottom <- apply(clu@assays$RNA@data[villus_bottom[villus_bottom %in% rownames(clu@assays$RNA@data)],], 2, median)
clu$villus_mid <- apply(clu@assays$RNA@data[villus_mid[villus_mid %in% rownames(clu@assays$RNA@data)],], 2, median)
clu$villus_top <- apply(clu@assays$RNA@data[villus_top[villus_top %in% rownames(clu@assays$RNA@data)],], 2, median)
DotPlot(clu, features = c('cvj', 'villus_bottom', 'villus_mid', 'villus_top'), col.min = 0, idents =c(grep('^RevSC', levels(clu), value = T) ,grep('^FCC', levels(clu), value = T),grep('^CC', levels(clu), value = T), 'EC-R11', 'EC-R29', 'EC-R5', 'EC-R17', 'EC-R2','EC-R20', 'EC-R18', 'EC-R6', "EC-R7"  ,"EC-R33", "EC-R15" , grep('^PC', levels(clu), value = T),grep('^TC', levels(clu), value = T),grep('^EE', levels(clu), value = T), grep('^GC', levels(clu), value = T),grep('^Fibro', levels(clu), value = T), grep('^Immune', levels(clu), value = T) )) & ggplot2::scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black')) & ggplot2::theme(text = element_text(size = 15),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15))
DotPlot(clu, features = c('cvj', 'villus_bottom', 'villus_mid', 'villus_top'), col.min = 0, idents =c(grep('^EC', levels(clu), value = T))) & ggplot2::scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black')) & ggplot2::theme(text = element_text(size = 15),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15))
DotPlot(clu, features = c('cvj', 'villus_bottom', 'villus_mid', 'villus_top'), col.min = 0, idents =levels(clu)[!levels(clu) %in% c(grep('Fibro', levels(clu), value = T), grep('Immune', levels(clu), value = T))]) & ggplot2::scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black')) & ggplot2::theme(text = element_text(size = 15),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15))
c(grep('^RevSC', levels(clu), value = T) ,grep('^FCC', levels(clu), value = T),grep('^CC', levels(clu), value = T), 'EC-W1', 'EC-W14', 'EC-W4' 'EC-W3', 'EC-W11', 'EC-W10', 'EC-W5', grep('^PC', levels(clu), value = T),grep('^TC', levels(clu), value = T),grep('^EE', levels(clu), value = T), grep('^GC', levels(clu), value = T),grep('^Fibro', levels(clu), value = T), grep('^Immune', levels(clu), value = T)
levels(clu) <- c(grep('^RevSC', levels(clu), value = T) ,grep('^FCC', levels(clu), value = T),grep('^CC', levels(clu), value = T), 'EC-W1', 'EC-W14', 'EC-W4' ,'EC-W3', 'EC-W11', 'EC-W10', 'EC-W5', grep('^PC', levels(clu), value = T),grep('^TC', levels(clu), value = T),grep('^EE', levels(clu), value = T), grep('^GC', levels(clu), value = T),grep('^Fibro', levels(clu), value = T), grep('^Immune', levels(clu), value = T))
levels(clu)[!levels(clu) %in% c(grep('^RevSC', levels(clu), value = T) ,grep('^FCC', levels(clu), value = T),grep('^CC', levels(clu), value = T), 'EC-W1', 'EC-W14', 'EC-W4' ,'EC-W3', 'EC-W11', 'EC-W10', 'EC-W5', grep('^PC', levels(clu), value = T),grep('^TC', levels(clu), value = T),grep('^EE', levels(clu), value = T), grep('^GC', levels(clu), value = T),grep('^Fibro', levels(clu), value = T), grep('^Immune', levels(clu), value = T))]
levels(clu)[!levels(clu) %in% c(grep('^RevSC', levels(clu), value = T) ,grep('^FCC', levels(clu), value = T),grep('^ICC', levels(clu), value = T), 'EC-W1', 'EC-W14', 'EC-W4' ,'EC-W3', 'EC-W11', 'EC-W10', 'EC-W5', 'EC-W2', 'EC-W7', grep('^PC', levels(clu), value = T),grep('^TC', levels(clu), value = T),grep('^EE', levels(clu), value = T), grep('^GC', levels(clu), value = T),grep('^Fibro', levels(clu), value = T), grep('^Immune', levels(clu), value = T))]
levels(clu) <- c(grep('^RevSC', levels(clu), value = T) ,grep('^FCC', levels(clu), value = T),grep('^ICC', levels(clu), value = T), 'EC-W1', 'EC-W14', 'EC-W4' ,'EC-W3', 'EC-W11', 'EC-W10', 'EC-W5', 'EC-W2', 'EC-W7', grep('^PC', levels(clu), value = T),grep('^TC', levels(clu), value = T),grep('^EE', levels(clu), value = T), grep('^GC', levels(clu), value = T),grep('^Fibro', levels(clu), value = T), grep('^Immune', levels(clu), value = T))
DotPlot(clu, features = c('cvj', 'villus_bottom', 'villus_mid', 'villus_top'), col.min = 0, idents =levels(clu)[!levels(clu) %in% c(grep('Fibro', levels(clu), value = T), grep('Immune', levels(clu), value = T))]) & ggplot2::scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black')) & ggplot2::theme(text = element_text(size = 15),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15))
DotPlot(clu, features = c('Alpi', 'Aoc1', 'Ccl25', 'Muc2', 'Agr2', 'Chga', 'Chgb', 'Dclk1', 'Trpm5','Defa17', 'Defa22', 'Lgr5', 'Olfm4', 'Clu', 'F3', 'Atg9b', 'Anxa1', 'Ly6a', 'Mki67','Pdgfra', 'Col1a1', 'Cd3g', 'Cd160', 'Lyz2', 'Csf1r'), col.min = 0) & scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black')) & theme(panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
DotPlot(clu, features = c('Alpi', 'Aoc1', 'Ccl25', 'Muc2', 'Agr2', 'Chga', 'Chgb', 'Dclk1', 'Trpm5','Defa17', 'Defa22', 'Lgr5', 'Olfm4', 'Clu', 'F3', 'Atg9b', 'Anxa1', 'Ly6a', 'Mki67','Pdgfra', 'Col1a1', 'Cd3g', 'Cd160', 'Lyz2', 'Csf1r'), col.min = 0) & scale_color_gradientn(colors = c('lightblue', 'yellow', 'red', 'black')) & theme(panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
clu$zonation2 <- colnames(clu@meta.data[,c('cvj', 'villus_bottom', 'villus_mid', 'villus_top')])[max.col(clu@meta.data[,c('cvj', 'villus_bottom', 'villus_mid', 'villus_top')], ties.method = 'first')]
entero <- clu[,WhichCells(clu, expression = collapse_idents == 'Enterocytes')]
Idents(clu) <- clu$wnn_clusters_idents
for(i in c('EC', 'GC', 'EE', 'PC', 'TC','ICC' ,'Immune', 'Fibro')){
clu <- SetIdent(clu, cells = WhichCells(clu, expression = wnn_clusters_idents %in% grep(paste('^', i, '.*', sep = ""), unique(clu$wnn_clusters_idents), value = T)), value = i)
}
clu <- RenameIdents(clu, 'ICC' = 'Crypt Cells', 'PC' = 'Paneth Cells', 'TC' = 'Tuft Cells', 'EE' = 'Enteroendocrine Cells', 'GC' = 'Goblet Cells', 'Fibro' = 'Fibroblasts', 'EC' = "Enterocytes")
levels(clu) <- c(grep('^RevSC.*', levels(clu), value = T) ,grep('^FCC.*', levels(clu), value = T),grep('^Crypt.*', levels(clu), value = T), grep('^Paneth.*', levels(clu), value = T),grep('^Tuft.*', levels(clu), value = T),grep('^Enteroendo.*', levels(clu), value = T), grep('^Goblet.*', levels(clu), value = T),grep('^Enterocytes.*', levels(clu), value = T), grep('^31.*', levels(clu), value = T),grep('^Fibro.*', levels(clu), value = T),grep('^Immune.*', levels(clu), value = T))
library(ggplot2)
clu$collapse_idents <- Idents(clu)
entero <- clu[,WhichCells(clu, expression = collapse_idents == 'Enterocytes')]
zonation_vs_time <- as.data.frame(table(entero$zonation2, entero$time))
colnames(zonation_vs_time) <- c('Identity', 'Time', 'Freq')
ggplot(zonation_vs_time, aes(fill = Time, y = Freq, x = Identity)) + geom_bar(position = 'dodge', stat = 'identity') + theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15)) + xlab('Identity')
ggplot(zonation_vs_time, aes(fill = Time, y = Freq, x = Identity)) + geom_bar(position = 'dodge', stat = 'identity') + theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15)) + xlab('Identity') + scale_fill_manual(values = c('#2166AC', '#B2182B'))
ggplot(zonation_vs_time, aes(fill = Time, y = Freq, x = Identity)) + geom_bar(position = 'stack', stat = 'identity') + theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15)) + xlab('Identity') + scale_fill_manual(values = c('#2166AC', '#B2182B'))
clucollap_vs_time <- as.data.frame(table(clu$collapse_idents, clu$time))
for(i in unique(clu$collapse_idents)){
clucollap_vs_time[which(clucollap_vs_time$Var1 == i),'prop'] <- clucollap_vs_time[which(clucollap_vs_time$Var1 == i),'Freq'] /length(WhichCells(clu, expression = collapse_idents == i))
}
colnames(clucollap_vs_time) <- c('Identity', 'Time', 'Freq', 'prop')
ggplot(clucollap_vs_time, aes(fill = Time, y = prop, x = Identity)) + geom_bar(position = 'stack', stat = 'identity') + theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15)) + xlab('Identity') + ylab('Temporal Proportion')
levels(clu)
clucollap_vs_time$Identity <- factor(clucollap_vs_time$Identity, levels = levels(clu))
ggplot(clucollap_vs_time, aes(fill = Time, y = prop, x = Identity)) + geom_bar(position = 'stack', stat = 'identity') + theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15)) + xlab('Identity') + ylab('Temporal Proportion')
ggplot(clucollap_vs_time, aes(fill = Time, y = prop, x = Identity)) + geom_bar(position = 'stack', stat = 'identity') + theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15)) + xlab('Identity') + ylab('Relative Temporal Proportion')
levels(clu)
levels(clu) <- c(grep('^RevSC.*', levels(clu), value = T) ,grep('^FCC.*', levels(clu), value = T),grep('^Crypt.*', levels(clu), value = T), grep('^Enterocytes.*', levels(clu), value = T),grep('^Paneth.*', levels(clu), value = T),grep('^Tuft.*', levels(clu), value = T),grep('^Enteroendo.*', levels(clu), value = T), grep('^Goblet.*', levels(clu), value = T), grep('^31.*', levels(clu), value = T),grep('^Fibro.*', levels(clu), value = T),grep('^Immune.*', levels(clu), value = T))
clucollap_vs_time$Identity <- factor(clucollap_vs_time$Identity, levels = levels(clu))
levels(clu) <- c(grep('^RevSC.*', levels(clu), value = T) ,grep('^FCC.*', levels(clu), value = T),grep('^Crypt.*', levels(clu), value = T), grep('^Enterocytes.*', levels(clu), value = T),grep('^Paneth.*', levels(clu), value = T),grep('^Tuft.*', levels(clu), value = T),grep('^Enteroendo.*', levels(clu), value = T), grep('^Goblet.*', levels(clu), value = T), grep('^31.*', levels(clu), value = T),grep('^Fibro.*', levels(clu), value = T),grep('^Immune.*', levels(clu), value = T))
levels(clu)
ggplot(clucollap_vs_time, aes(fill = Time, y = prop, x = Identity)) + geom_bar(position = 'stack', stat = 'identity') + theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15)) + xlab('Identity') + ylab('Relative Temporal Proportion')
Idents(clu) <- clu$new_wnn_clusters
levels(clu) <- c(grep('^RevSC', levels(clu), value = T) ,grep('^FCC', levels(clu), value = T),grep('^ICC', levels(clu), value = T), 'EC-W1', 'EC-W14', 'EC-W4' ,'EC-W3', 'EC-W11', 'EC-W10', 'EC-W5', 'EC-W2', 'EC-W7', grep('^PC', levels(clu), value = T),grep('^TC', levels(clu), value = T),grep('^EE', levels(clu), value = T), grep('^GC', levels(clu), value = T),grep('^Fibro', levels(clu), value = T), grep('^Immune', levels(clu), value = T))
cluwnn_vs_time <- as.data.frame(table(clu$new_wnn_clusters, clu$time))
for(i in unique(clu$new_wnn_clusters)){
cluwnn_vs_time[which(cluwnn_vs_time$Var1 == i),'prop'] <- cluwnn_vs_time[which(cluwnn_vs_time$Var1 == i),'Freq'] /length(WhichCells(clu, expression = new_wnn_clusters == i))
}
cluwnn_vs_time$Identity <- factor(cluwnn_vs_time$Identity, levels = levels(clu))
levels(clu)
View(cluwnn_vs_time)
colnames(cluwnn_vs_time) <- c('Identity', 'Time', 'Freq', 'prop')
cluwnn_vs_time$Identity <- factor(cluwnn_vs_time$Identity, levels = levels(clu))
ggplot(cluwnn_vs_time, aes(fill = Time, y = prop, x = Identity)) + geom_bar(position = 'stack', stat = 'identity') + theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15)) + xlab('Identity') + ylab('Relative Temporal Proportion')
ggplot(cluwnn_vs_time, aes(fill = Time, y = prop, x = Identity)) + geom_bar(position = 'stack', stat = 'identity') + theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15)) + xlab('Identity') + ylab('Relative Temporal Proportion') + scale_fill_manual(values = c('#2166ac', '#b2182b'))
Idents(clu) <- clu$collapse_idents
levels(clu) <- c(grep('^RevSC.*', levels(clu), value = T) ,grep('^FCC.*', levels(clu), value = T),grep('^Crypt.*', levels(clu), value = T), grep('^Enterocytes.*', levels(clu), value = T),grep('^Paneth.*', levels(clu), value = T),grep('^Tuft.*', levels(clu), value = T),grep('^Enteroendo.*', levels(clu), value = T), grep('^Goblet.*', levels(clu), value = T), grep('^31.*', levels(clu), value = T),grep('^Fibro.*', levels(clu), value = T),grep('^Immune.*', levels(clu), value = T))
clucollap_vs_time <- as.data.frame(table(clu$collapse_idents, clu$time))
for(i in unique(clu$collapse_idents)){
clucollap_vs_time[which(clucollap_vs_time$Var1 == i),'prop'] <- clucollap_vs_time[which(clucollap_vs_time$Var1 == i),'Freq'] /length(WhichCells(clu, expression = collapse_idents == i))
}
colnames(clucollap_vs_time) <- c('Identity', 'Time', 'Freq', 'prop')
ggplot(clucollap_vs_time, aes(fill = Time, y = prop, x = Identity)) + geom_bar(position = 'stack', stat = 'identity') + theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15)) + xlab('Identity') + ylab('Temporal Proportion')
clucollap_vs_time$Identity <- factor(clucollap_vs_time$Identity, levels = levels(clu))
for(i in unique(clu$collapse_idents)){
clucollap_vs_time[which(clucollap_vs_time$Var1 == i),'prop'] <- clucollap_vs_time[which(clucollap_vs_time$Var1 == i),'Freq'] /length(WhichCells(clu, expression = collapse_idents == i))
}
colnames(clucollap_vs_time) <- c('Identity', 'Time', 'Freq', 'prop')
ggplot(clucollap_vs_time, aes(fill = Time, y = prop, x = Identity)) + geom_bar(position = 'stack', stat = 'identity') + theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15)) + xlab('Identity') + ylab('Temporal Proportion')
clucollap_vs_time <- as.data.frame(table(clu$collapse_idents, clu$time))
for(i in unique(clu$collapse_idents)){
clucollap_vs_time[which(clucollap_vs_time$Var1 == i),'prop'] <- clucollap_vs_time[which(clucollap_vs_time$Var1 == i),'Freq'] /length(WhichCells(clu, expression = collapse_idents == i))
}
colnames(clucollap_vs_time) <- c('Identity', 'Time', 'Freq', 'prop')
clucollap_vs_time$Identity <- factor(clucollap_vs_time$Identity, levels = levels(clu))
ggplot(clucollap_vs_time, aes(fill = Time, y = prop, x = Identity)) + geom_bar(position = 'stack', stat = 'identity') + theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15)) + xlab('Identity') + ylab('Temporal Proportion')
ggplot(clucollap_vs_time, aes(fill = Time, y = prop, x = Identity)) + geom_bar(position = 'stack', stat = 'identity') + theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15)) + xlab('Identity') + ylab('Temporal Proportion') + scale_fill_manual(values = c('#2166ac', '#b2182b'))
df <- data.frame(matrix(nrow = length(levels(clu)), ncol = 2))
rownames(df) <- levels(clu)
colnames(df) <- c('RNA.weight', 'peaks.weight')
for(cluster in levels(clu)){
df[cluster,'RNA.weight'] <- median(clu@meta.data[WhichCells(clu, idents =cluster),'RNA.weight'])
df[cluster,'peaks.weight'] <- median(clu@meta.data[WhichCells(clu, idents =cluster),'peaks.weight'])
}
df$clusters <- rownames(df)
library(reshape2)
dfm <- melt(df[,c('RNA.weight', 'peaks.weight', 'clusters')], id.vars = 3)
colnames(dfm) <- c('clusters', 'Modality', 'value')
dfm$clusters <- factor(dfm$clusters, levels = levels(clu))
ggplot(dfm, aes(x = clusters, y = value, fill = Modality))+ geom_bar(position = 'dodge', stat = 'identity') + theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15)) + xlab('Identity')
Idents(clu) <- clu$new_wnn_clusters
levels(clu) <- c(grep('^RevSC', levels(clu), value = T) ,grep('^FCC', levels(clu), value = T),grep('^ICC', levels(clu), value = T), 'EC-W1', 'EC-W14', 'EC-W4' ,'EC-W3', 'EC-W11', 'EC-W10', 'EC-W5', 'EC-W2', 'EC-W7', grep('^PC', levels(clu), value = T),grep('^TC', levels(clu), value = T),grep('^EE', levels(clu), value = T), grep('^GC', levels(clu), value = T),grep('^Fibro', levels(clu), value = T), grep('^Immune', levels(clu), value = T))
df <- data.frame(matrix(nrow = length(levels(clu)), ncol = 2))
rownames(df) <- levels(clu)
colnames(df) <- c('RNA.weight', 'peaks.weight')
for(cluster in levels(clu)){
df[cluster,'RNA.weight'] <- median(clu@meta.data[WhichCells(clu, idents =cluster),'RNA.weight'])
df[cluster,'peaks.weight'] <- median(clu@meta.data[WhichCells(clu, idents =cluster),'peaks.weight'])
}
df$clusters <- rownames(df)
library(reshape2)
dfm <- melt(df[,c('RNA.weight', 'peaks.weight', 'clusters')], id.vars = 3)
colnames(dfm) <- c('clusters', 'Modality', 'value')
dfm$clusters <- factor(dfm$clusters, levels = levels(clu))
ggplot(dfm, aes(x = clusters, y = value, fill = Modality))+ geom_bar(position = 'dodge', stat = 'identity') + theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15)) + xlab('Identity')
ggplot(dfm, aes(x = clusters, y = value, fill = Modality))+ geom_bar(position = 'dodge', stat = 'identity') + theme(text = element_text(size = 15, family = 'Helvetica'),axis.line = element_line(size = 1),axis.ticks = element_line(size = 1),panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15)) + xlab('Identity') + scale_fill_manual(values = c('#2166ac', '#b2182b'))
save(clu, file = 'cluall_m5_0dpi_2dpi_dblt_filtd_rna_clustered_idented_peaks_recalled_genesLinked_chromvar_weighted_macs2peaks_embed_clustered_idented_cluster_formatted.RData')
