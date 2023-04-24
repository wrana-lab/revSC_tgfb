setwd("/media/shyam/external/multiome_clu_2/multi5")
library(Seurat)
library(Signac)
setwd("/media/shyam/external/multiome")
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(EnsDb.Mmusculus.v79)
library(Signac)
for(i in c('C1_A9','C3_A11')){
val <- paste(i,"/outs/atac_peaks.bed", sep = "")
peaks <- read.table(file = val, col.names = c('chr', 'start', 'end'))
peaks <- makeGRangesFromDataFrame(peaks)
assign(paste0("peaks_", i), peaks)
}
combined.peaks <- reduce(x = c(peaks_C1_A9, peaks_C3_A11))
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
for(i in c('C1_A9','C3_A11')){
val <- paste(i,"/outs/filtered_feature_bc_matrix/", sep = "")
print(val)
assign(paste0("obj_", i), Read10X(data.dir=val)$`Gene Expression`)
}
rm(i)
rm(val)
for(i in grep("obj_", objects(), value = T)){
print(i)
assign(paste0("", i), CreateSeuratObject(counts = get(i)))
}
for(i in grep("obj_", objects(), value = T)){
print(i)
objj <- get(i)
scrubs <- read.csv(file = paste("scrublets/",stringr::str_to_lower(sub("obj_", "", i)), '_scrublets.csv', sep = ""), sep = ',', header = FALSE)
rownames(scrubs) <- colnames(objj)
objj$scrub <- scrubs[,1]
assign(paste0("", i), objj)
}
#for(i in grep("obj_", objects(), value = T)){
#  objj <- get(i)
#  objj <- subset(objj, subset = nFeature_spliced > 200)
#  assign(paste0("", i), objj)
#}
library(stringr)
setwd("/media/shyam/external/multiome_clu_2/multi5")
blacklist_mm39 <- read.csv(file = 'mm39_blacklist_frommm10_liftover.csv', header = F)
blacklist_mm39 <- as.data.frame(blacklist_mm39)
colnames(blacklist_mm39)
blacklist_mm39 <- makeGRangesFromDataFrame(blacklist_mm39, start.field = 'V2', end.field = 'V3', seqnames.field = 'V1', starts.in.df.are.0based = T)
blacklist_mm39 <- as.data.frame(blacklist_mm39)
blacklist_mm39$seqnames <-   str_replace(blacklist_mm39$seqnames, 'chr', "")
blacklist_mm39 <- makeGRangesFromDataFrame(blacklist_mm39, keep.extra.columns = T)
setwd("/media/shyam/external/GRCm39/genes")
gtf <- rtracklayer::import('genes.gtf')
setwd("/media/shyam/external/multiome")
gene.coords <- gtf
gene.coords$gene_id
annotations <- gene.coords
for(i in c('C1_A9','C3_A11')){
val <- paste(i,"/outs/per_barcode_metrics.csv", sep = "")
print(val)
metadata <- read.csv(file = val, header = T, row.names = 1)
print(i)
objj <- get(paste0('obj_', i, sep = ''))
#objj <- CreateSeuratObject(counts = rna_counts, meta.data = metadata)
objj <- AddMetaData(objj, metadata)
frag.file <- paste(i,"/outs/atac_fragments.tsv.gz", sep = "")
frags <- CreateFragmentObject(path = frag.file)
fragcounts <- FeatureMatrix(fragments = frags, features = combined.peaks, cells = objj$gex_barcode)
#annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
#seqlevelsStyle(annotations) <- 'UCSC'
#genome(annotations) <- "mm10"
fragatc <- CreateChromatinAssay(
counts = fragcounts,
sep = c(":", "-"),
genome = 'mm39',
fragments = frags,
annotation = annotations,
)
#fragatac <- CreateChromatinAssay(fragcounts, fragments = frags)
fragatc <- RenameCells(fragatc, new.names = colnames(objj))
objj[["ATAC"]] <- fragatc
DefaultAssay(objj) <- 'ATAC'
objj <- NucleosomeSignal(object = objj)
objj <- TSSEnrichment(object = objj, fast = FALSE)
objj$pct_reads_in_peaks <- objj$atac_peak_region_fragments / objj$atac_fragments
objj$blacklist_fraction <- FractionCountsInRegion(objj, assay = 'ATAC', regions = blacklist_mm39)
objj$blacklist_ratio <- objj$blacklist_fraction / objj$atac_peak_region_fragments
#mito_genes <- rownames(x = objj@assays$RNA@counts)[grep(pattern = "^mt-", x = rownames(x =objj@assays$RNA@counts)) ]
#percent.mito <- colSums(objj@assays$RNA@counts[mito_genes, ])/colSums(objj@assays$RNA@counts)
#objj$percent_mito <- percent.mito
DefaultAssay(objj) <- 'RNA'
assign(paste0("obj_", i), objj)
}
df <- data.frame(matrix(ncol=5, nrow = length(grep("obj_", objects(), value = T))))
#x <- c("obj_name", 'initial_cells', "after_feature_filter", "after_lib_size", "after_mito_filter")
x <- c("obj_name", 'initial_cells', "after_mito_filter", "after_feature_filter", "after_lib_size")
colnames(df) <- x
df$obj_name <- grep("obj_", objects(), value = T)
for(i in grep("obj_", objects(), value = T)){
print(i)
objj <- get(i)
df$initial_cells[df$obj_name == i] = ncol(objj)
mito_genes <- rownames(x = objj@assays$RNA@counts)[grep(pattern = "^mt-", x = rownames(x =objj@assays$RNA@counts)) ]
percent.mito <- colSums(objj@assays$RNA@counts[mito_genes, ])/colSums(objj@assays$RNA@counts)
objj$percent_mito <- percent.mito
assign(paste0("", i), objj)
}
head(obj_C1_A9$percent_mito)
library(scater)
library(scran)
library(scuttle)
library(scDblFinder)
qcer <- function(a, obj, query, q, b){
f0 <- qplot(obj[[query]][,1], geom = 'histogram')
filt <- as.data.frame(seq(0,a, by = b))
filt['filt'] = filt[,1]
filt['val'] = 0
for(i in 1:nrow(filt)){
filt[i,'val'] <- length(obj[[query]][obj[[query]] > filt[i,'filt']])
}
filt['diff'] = 0
for(i in 2:nrow(filt)){
filt[i,'diff'] <- filt[i,'val'] -  filt[i-1,'val']
}
print(filt)
i = deparse(substitute(obj))
q[paste(i, query, 'filter_dataframe', sep='_')] <- list(filt)
f1 = ggplot(filt, aes(x = filt, y= val)) + geom_point()+ geom_line()
f2 = ggplot(filt, aes(x = filt, y= diff)) + geom_point()+ geom_line()
q[paste(i, query, 'filter_val', sep='_')] <-list(f1)
q[paste(i, query, 'filter_diff', sep='_')] <-list(f2)
q[paste(i, query, 'qplot', sep='_')] <-list(f0)
r <- list()
r <- append(r, q)
return(r)
}
qc_obj <- function(objj, featurefilt, countfilt, nucfilt, tssfilt, peaksfilt, mitofilt, frips, blacklistlimit, r){
i = deparse(substitute(objj))
objj <- objj[,rownames(na.omit(objj@meta.data))]
objsce <- as.SingleCellExperiment(objj)
keep_feature <- rowSums(counts(objsce) > 0) >= 3
rna <- GetAssay(objj, assay = 'RNA')
rna <- rna[keep_feature,]
rna <- CreateAssayObject(rna)
objj[['RNA']] <- rna
DefaultAssay(objj) <- 'RNA'
objsce <- as.SingleCellExperiment(objj)
#35%
objj$ncountfeat <- ifelse(objj$nFeature_RNA<featurefilt, T, F)
objj$ncountfilt <- ifelse(objj$nCount_RNA<countfilt, T, F)
#filter on high, means more mononuc than nuc free frags
objj$nuc_filt <- ifelse(objj$nucleosome_signal<nucfilt, T, F)
#filter on low
objj$tss_filt <- ifelse(objj$TSS.enrichment<tssfilt, T, F)
#filter on low
objj$peaks_filt <- ifelse(objj$atac_peak_region_fragments<peaksfilt, T, F)
#filter on low
#objj$reads_filt <- ifelse(objj$pct_reads_in_peaks<quantile(objj$pct_reads_in_peaks, probs = seq(0.05, 1.0, by = .05))[4], T, F)
#filter on high
#objj$blacklist_filt <- ifelse(objj$blacklist_ratio<quantile(objj$blacklist_ratio, probs = seq(0.05, 1.0, by = .05))[16], T, F)
#keeped = WhichCells(objj, expression = ncountfeat == F & ncountfilt ==F & percent_mito <mitofilt & nuc_filt == T & tss_filt == F & peaks_filt == F & pct_reads_in_peaks > frips & blacklist_ratio < blacklistlimit)
keeped = WhichCells(objj, expression = nFeature_RNA>featurefilt & nCount_RNA>countfilt & percent_mito <mitofilt & nucleosome_signal<nucfilt & TSS.enrichment>tssfilt & atac_peak_region_fragments>peaksfilt & pct_reads_in_peaks > frips & blacklist_ratio < blacklistlimit)
#keeped = WhichCells(objj, expression = ncountfeat == F & ncountfilt ==F & percent_mito <30)
objj$keeped <- is.na(match(colnames(objj), keeped))
swag <- FeatureScatter(objj, feature1 = 'nCount_RNA', feature2='percent_mito', group.by = 'keeped')
r[paste(i, 'filtered_mito_vs_ncount', sep = '_')] <- list(swag)
swag <- FeatureScatter(objj, feature1 = 'nCount_RNA', feature2='nFeature_RNA', group.by = 'keeped')
r[paste(i, 'filtered_feature_vs_ncount', sep = '_')] <- list(swag)
swag <- FeatureScatter(objj, feature1 = 'nFeature_RNA', feature2='percent_mito', group.by = 'keeped')
r[paste(i, 'filtered_mito_vs_feature', sep = '_')] <- list(swag)
swag <- VlnPlot(objj, features = c('nucleosome_signal', 'TSS.enrichment','pct_reads_in_peaks', 'blacklist_ratio' ), group.by = 'keeped')
r[paste(i, 'atac_metrics', sep = '_')] <- list(swag)
objsce <- scDblFinder(objsce)
objj$Doublet <- objsce$scDblFinder.class
singlets = WhichCells(objj, expression = Doublet == 'singlet')
keeped = keeped[keeped %in% singlets]
objsce <- objsce[,keeped]
obj_clus <- quickCluster(objsce)
objsce <- computeSumFactors(objsce, clusters = obj_clus)
objsce <- logNormCounts(objsce)
objj_2 <- as.Seurat(objsce)
objj <- objj[,keeped]
objj[['RNA']] <- objj_2[['RNA']]
swag <- FeatureScatter(objj, feature1 = 'nCount_RNA', feature2='nFeature_RNA', group.by = 'Doublet')
r[paste(i, 'doublets_nFeature_vs_count', sep = '_')] <- list(swag)
objj2 <- subset(objj, subset = Doublet=="singlet")
objlist <- list()
objlist['qclist'] <- list(r)
objlist['obj'] <- list(objj)
return(objlist)
}
mitofilt <- 0.50
frips <- 0.15
blacklistlimit <- 0.05
mitofilt <- 0.50
frips <- 0.15
blacklistlimit <- 0.05
q = list()
r <- list()
b <- 250
qplot(obj_C1_A9$nFeature_RNA, geom = 'histogram')
a = 4000
q = qcer(a, obj_C1_A9, 'nFeature_RNA', q, b)
q$obj_C1_A9_nFeature_RNA_filter_diff
################################
featurefilt <- 250
#########################3
qplot(obj_C1_A9$nCount_RNA, geom = 'histogram')
a = 10000
q = qcer(a, obj_C1_A9, 'nCount_RNA', q, b)
q$obj_C1_A9_nCount_RNA_filter_diff
################################
countfilt <- 250
################################
qplot(obj_C1_A9$nucleosome_signal, geom = 'histogram')
a = 1
b = 0.025
q = qcer(a, obj_C1_A9, 'nucleosome_signal', q, b)
################################
nucfilt <- 4
################################
qplot(obj_C1_A9$TSS.enrichment, geom = 'histogram')
a = 20
b = 0.25
q = qcer(a, obj_C1_A9, 'TSS.enrichment', q, b)
q$obj_C1_A9_TSS.enrichment_filter_diff
################################
tssfilt <- 2.25
################################
qplot(obj_C1_A9$atac_peak_region_fragments, geom = 'histogram')
a = 20000
b = 250
q = qcer(a, obj_C1_A9, 'atac_peak_region_fragments', q, b)
q$obj_C1_A9_atac_peak_region_fragments_filter_diff
################################
peaksfilt <- 100
################################
obs <- qc_obj(obj_C1_A9, featurefilt, countfilt, nucfilt, tssfilt, peaksfilt, mitofilt, frips, blacklistlimit, r)
obj_C1_A9 <- obs$obj
r <- obs$qclist
b <- 250
qplot(obj_C3_A11$nFeature_RNA, geom = 'histogram')
a = 4000
q = qcer(a, obj_C3_A11, 'nFeature_RNA', q, b)
q$obj_C3_A11_nFeature_RNA_filter_diff
################################
featurefilt <- 250
################################
qplot(obj_C3_A11$nCount_RNA, geom = 'histogram')
a = 10000
q = qcer(a, obj_C3_A11, 'nCount_RNA', q, b)
q$obj_C3_A11_nCount_RNA_filter_diff
################################
countfilt <- 250
################################
qplot(obj_C3_A11$nucleosome_signal, geom = 'histogram')
a = 1
b = 0.025
q = qcer(a, obj_C3_A11, 'nucleosome_signal', q, b)
################################
nucfilt <- 4
################################
qplot(obj_C3_A11$TSS.enrichment, geom = 'histogram')
a = 20
b = 0.25
q = qcer(a, obj_C3_A11, 'TSS.enrichment', q, b)
q$obj_C3_A11_TSS.enrichment_filter_diff
################################
tssfilt <- 2.00
################################
qplot(obj_C3_A11$atac_peak_region_fragments, geom = 'histogram')
a = 20000
b = 250
q = qcer(a, obj_C3_A11, 'atac_peak_region_fragments', q, b)
q$obj_C3_A11_atac_peak_region_fragments_filter_diff
################################
peaksfilt <- 100
################################
obs <- qc_obj(obj_C3_A11, featurefilt, countfilt, nucfilt, tssfilt, peaksfilt, mitofilt, frips, blacklistlimit, r)
obj_C3_A11 <- obs$obj
r <- obs$qclist
df$afterfilter <- NA
for(i in grep("obj_", objects(), value = T)){
objj <- get(i)
df$afterfilter[df$obj_name == i] = ncol(objj)
}
df$total_filtered <- df$initial_cells - df$afterfilter
#masterlist <- mget(ls(pattern="obj_"))
rm(obj_clus)
for(i in grep("obj_", objects(), value = T)){
objj <- get(i)
objj$batch = i
assign(paste0("", i), objj)
}
for(i in grep("obj_", objects(), value = T)){
objj <- get(i)
DefaultAssay(objj) <- 'ATAC'
objj <- RunTFIDF(objj, assay = 'ATAC')
objj <- FindTopFeatures(objj, min.cutoff = 10, assay = 'ATAC')
objj <- RunSVD(objj, assay = 'ATAC')
DefaultAssay(objj) <- 'RNA'
assign(paste0("", i), objj)
}
clu <- merge(x = obj_C3_A11, y = obj_C1_A9, merge.data = T)
#clu$time <- ifelse(clu$batch == 'obj_C1_A9' | clu$batch == 'obj_c4_hnn' | clu$batch == 'obj_c5_hnn' | clu$batch == 'obj_c4_hl5' | clu$batch == 'obj_c5_hl5', '0dpi', ifelse(clu$batch == 'obj_C2_A10', '1dpi', ifelse(clu$batch == 'obj_C3_A11', '2dpi', ifelse(clu$batch == 'obj_C4_A12', '4dpi', '3dpi')))
clu$time <- ifelse(clu$batch == 'obj_C1_A9', '0dpi', ifelse(clu$batch == 'obj_C2_A10', '1dpi', ifelse(clu$batch == 'obj_C3_A11', '2dpi', ifelse(clu$batch == 'obj_C4_A12', '4dpi', '3dpi'))))
setwd("/media/shyam/external/multiome_clu_2/multi5")
save(clu, file = 'cluall_multi5_0dpi_2dpi_preprocessed_dblt_filtd.RData')
save(r, file = 'qc_clu_multi5_0dpi_2dpi_dblt_filtd.RData')
save(df, file = 'df_clu_multi5_0dpi_2dpi_dblt_filtd.RData')
save(q, file = 'preqc_clu_multi5_0dpi_2dpi_dblt_filtd.RData')
savehistory("/media/shyam/external/multiome_clu_2/multi5/cluall_multi5_0dpi_2dpi_preprocessed_dblt_filtd_real.R")
