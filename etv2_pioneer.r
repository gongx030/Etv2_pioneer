
# [2017-11-20] Prepare the ATAC-seq data
project <- 'Etv2_MEFs'
project.dir <- sprintf('%s/projects/%s', Sys.getenv('SHARED'), project)
bam.files <- c(
	'D2_EB_DOX' = sprintf('%s/mESCs_bams/d3_EB_Dox_1_S3.rmdup.bam', project.dir),
	'D3_EB_DOX' = sprintf('%s/mESCs_bams/D3EB_F_P_6h_Dox-1_S12.rmdup.bam',project.dir),
	'D3_EB' 		= sprintf('%s/mESCs_bams/D3EB_F_P_no_Dox-1_S10.rmdup.bam',project.dir),
	'D2_EB' 		= sprintf('%s/mESCs_bams/d3_EB_no_dox_1_S1.rmdup.bam', project.dir),
	'ESC_DOX' 	= sprintf('%s/mESCs_bams/ES_6h_Dox-1_S8.rmdup.bam', project.dir),
	'ESC' 			= sprintf('%s/mESCs_bams/ES_no_Dox-1_S6.rmdup.bam', project.dir),
	'D2_EB_DOX' = sprintf('%s/mESCs_bams/d3_EB_Dox_2_S4.rmdup.bam', project.dir),
	'D3_EB_DOX' = sprintf('%s/mESCs_bams/D3EB_F_P_6h_Dox-2_S13.rmdup.bam', project.dir),
	'D3_EB' 		= sprintf('%s/mESCs_bams/D3EB_F_P_no_Dox-2_S11.rmdup.bam', project.dir),
	'D2_EB' 		= sprintf('%s/mESCs_bams/d3_EB_no_dox_2_S2.rmdup.bam', project.dir),
	'ESC_DOX' 	= sprintf('%s/mESCs_bams/ES_6h_Dox-2_S9.rmdup.bam', project.dir),
	'ESC' 			= sprintf('%s/mESCs_bams/ES_no_Dox-2_S7.rmdup.bam', project.dir),
	'MEF_DOX' 	= sprintf('%s/lab_and_public_MEFs_bams/3hDox1.rmdup.bam', project.dir),
	'MEF_DOX' 	= sprintf('%s/lab_and_public_MEFs_bams/3hDox2.rmdup.bam', project.dir),
	'MEF' 			= sprintf('%s/lab_and_public_MEFs_bams/noDox1.rmdup.bam', project.dir),
	'MEF' 			= sprintf('%s/lab_and_public_MEFs_bams/noDox2.rmdup.bam', project.dir),
	'MEF_SRR1930168' = sprintf('%s/lab_and_public_MEFs_bams/SRR1930168.rmdup.bam', project.dir),
	'mesoderm' = '/panfs/roc/scratch/gongx030/ncbi/sra/ENCFF/ENCFF31/ENCFF315DJU/ENCFF315DJU.bam',	# ENCODE meosderm E11.5
	'mesoderm' = '/panfs/roc/scratch/gongx030/ncbi/sra/ENCFF/ENCFF59/ENCFF590IWJ/ENCFF590IWJ.bam', 	# ENCODE mesoderm E11.5
	'Tie2pos' = '/panfs/roc/scratch/gongx030/ncbi/sra/ERR10/ERR1079/ERR1079314/ERR1079314.dedup.bam' # Tie2+ from E12.5
)
d <- data.frame(bam.file = bam.files, group = names(bam.files))
d <- transform(d, coverage.file = gsub('.bam$', '.bed', bam.file), bam.index.file = sprintf('%s.bai', bam.file))
table(file.exists(d[, 'bam.file'])) # make sure all files exist

# [2017-11-20] this is a union DHS sites of all ENCODE mouse ATAC-seq/DNase-seq sites
dhs.bed.file <- sprintf('%s/projects/sc_atac/mouse_dhs.300.bed', Sys.getenv('SHARED'))

# [2017-11-26] index the BAM files, if the index is not available; bedtools multicov only works for indexed BAM files
table(file.exists(d[, 'bam.index.file']))
library(parallel); mclapply(which(!file.exists(d[, 'bam.index.file'])), function(i){
	command <- sprintf('samtools index %s', d[i, 'bam.file'])
	cat(sprintf('[%s] %s\n', Sys.time(), command)); system(command)
}, mc.cores = 2)

# [2017-11-20] generate depth coverage for each BAM file, based on the input peaks
# run this on lab queues, since it may take >20 mins for large BAM files
table(file.exists(d[, 'coverage.file']))
library(parallel); mclapply(which(!file.exists(d[, 'coverage.file'])), function(i){
	bam.index.file <- sprintf('%s.bai', d[i, 'bam.file'])
	command <- sprintf('bedtools multicov -bams %s -bed %s > %s', d[i, 'bam.file'], dhs.bed.file, d[i, 'coverage.file'])
	cat(sprintf('[%s] %s\n', Sys.time(), command)); system(command)
}, mc.cores = 3)


# [2017-11-20] read the read count matrix (version=20171120a)
# [2017-11-26] Adding one ENCODE E11.5 mesoderm samples (ENCFF590IWJ), and endothelial cell (Tie2+ from E12.5, ERR1079314) (version=20171126a)
# [2017-11-26] Adding another ENCODE E11.5 mesoderm samples (ENCFF315DJU) (version=20171126b)

#X.file <- sprintf('%s/readcount_version=20171120a.RData', project.dir); peaks.file <- sprintf('%s/peaks_version=20171120a.RData', project.dir); motif_ix.file <- sprintf('%s/motif_ix.mouse_dhs.300.RData', project.dir)
#X.file <- sprintf('%s/readcount_version=20171126a.RData', project.dir); peaks.file <- sprintf('%s/peaks_version=20171126a.RData', project.dir); motif_ix.file <- sprintf('%s/motif_ix.mouse_dhs.300_version=20171126a.RData', project.dir)
X.file <- sprintf('%s/readcount_version=20171126b.RData', project.dir); peaks.file <- sprintf('%s/peaks_version=20171126b.RData', project.dir); motif_ix.file <- sprintf('%s/motif_ix.mouse_dhs.300_version=20171126b.RData', project.dir)

X <- as.matrix(do.call('cbind', lapply(d[, 'coverage.file'], function(file){
	cat(sprintf('[%s] reading %s\n', Sys.time(), file))
	read.table(file, header = FALSE, sep = '\t')[, 4]
})))	# takes a few mins to read; save the intermediate results
peaks <- read.table(dhs.bed.file, header = FALSE, sep = '\t')
rownames(X) <- rownames(peaks) <- sprintf('%s:%d-%d', peaks[, 1], peaks[, 2], peaks[, 3])
class(X) <- 'integer'
i <- rowSums(X > 1) > 0	# only include the peaks that are present in at least one sample
X <- X[i, ]; peaks <- peaks[i, ]
library(GenomicRanges); peaks <- GRanges(seqnames = peaks[, 1], range = IRanges(peaks[, 2], peaks[, 3]))
save(X, file = X.file)
save(peaks, file = peaks.file)

# [2017-11-20] chromVAR analysis on combined ATAC-seq data
X <- local({get(load(X.file))})
peaks <- local({get(load(peaks.file))})
library(chromVARmotifs) # https://github.com/GreenleafLab/chromVARmotifs
library(chromVAR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)
library(BiocParallel)
library(SummarizedExperiment)
register(MulticoreParam(4)) # Use 8 cores
motif.set <- 'homer_pwms'
#motif.set <- 'mouse_pwms_v2'
data(list = motif.set, package = 'chromVARmotifs')
rse <- SummarizedExperiment(assays = SimpleList(counts = X), rowRanges = peaks, colData = d)
rse <- addGCBias(rse, genome = BSgenome.Mmusculus.UCSC.mm10)
#motif_ix <- matchMotifs(get(motif.set), peaks, genome = 'mm10')
#save(motif_ix, file = motif_ix.file)
motif_ix <- local({get(load(motif_ix.file))})
dev <- computeDeviations(object = rse, annotations = motif_ix)
v <- computeVariability(dev)
												                                  
# visualize the TF variance by heatmap
#m <- !d[, 'group'] %in% 'mesoderm'
m <- rep(TRUE, nrow(d))
Z <- assay(dev)[order(v[, 'p_value_adj'])[1:130], m]	# choose the top most variable TFs
colnames(Z) <- d[m, 'group']
library(pheatmap); library(gplots)
col <- colorpanel(100, low = 'blue', mid = 'white', high = 'red')
breaks <- c(-10, seq(-0.25, 0.25, length.out = 99), 10)
pheatmap(Z, cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = TRUE, breaks = breaks, color = col, cellheight = 7, cellwidth = 30, fontsize_col = 15, fontsize_row = 7, show_colnames = TRUE, border_color = 'black', annotation_names_row = FALSE, annotation_names_col = FALSE)


# visualize the sample PCA plot
library(irlba); y <- irlba(assay(dev)[, m], nu = 1, nv = 2)$v
col <- c(
	'D2_EB' = 'green', 			'D2_EB_DOX' = 'green', 
	'D3_EB' = 'darkgreen', 	'D3_EB_DOX' = 'darkgreen', 	
	'ESC' = 'yellow', 			'ESC_DOX' = 'yellow', 
	'MEF' = 'blue', 				'MEF_DOX' = 'blue', 	'MEF_SRR1930168' = 'blue',
	'Tie2pos' = 'black', 			'mesoderm' = 'lightgreen'
)
pch <- c(
	'D2_EB' = 21,	'D2_EB_DOX' = 22, 			
	'D3_EB' = 21,	'D3_EB_DOX' = 22,
	'ESC' = 21, 'ESC_DOX' = 22, 
	'MEF' = 21, 'MEF_DOX' = 22, 'MEF_SRR1930168' = 21,
	'Tie2pos' = 21, 			'mesoderm' = 21
)
plot(y[, 1], y[, 2], col = 'black', bg = col[d[m, 'group']], pch = pch[d[m, 'group']], cex = 2)


x1 <- rowMeans(assay(dev)[, d[, 'group'] == 'MEF'])
x2 <- rowMeans(assay(dev)[, d[, 'group'] == 'MEF_DOX'])
yy <- sort(x2 - x1)
i <- yy > 0.04 
par(mar = c(10, 10, 10, 10)); plot(1:length(yy), yy, pch = 21, bg = 'gray', col = 'black', cex = 1, xaxt = 'n', ylim = c(min(yy), max(yy) * 1.5), yaxt = 'n', xlab = '', ylab = '')
library(wordcloud); textplot(which(i), yy[i], gsub('(.+?)\\/.+', '\\1', names(yy)[i]), xpd = TRUE, new = FALSE, adj = 1, cex = 1.25)


x1 <- rowMeans(assay(dev)[, d[, 'group'] == 'D2_EB'])
x2 <- rowMeans(assay(dev)[, d[, 'group'] == 'D2_EB_DOX'])
yy <- sort(x2 - x1)
i <- yy > 0.02 
par(mar = c(10, 10, 10, 10)); plot(1:length(yy), yy, pch = 21, bg = 'gray', col = 'black', cex = 1, xaxt = 'n', ylim = c(min(yy), max(yy) * 1.5), yaxt = 'n', xlab = '', ylab = '')
library(wordcloud); textplot(which(i), yy[i], gsub('(.+?)\\/.+', '\\1', names(yy)[i]), xpd = TRUE, new = FALSE, adj = 1, cex = 1.25)

z1 <- rowMeans(assay(dev)[, d[, 'group'] == 'MEF_DOX']) - rowMeans(assay(dev)[, d[, 'group'] == 'MEF']) 
z2 <- rowMeans(assay(dev)[, d[, 'group'] == 'D2_EB_DOX']) - rowMeans(assay(dev)[, d[, 'group'] == 'D2_EB']) 
bg <- rep('gray', length(z1)); cex <- rep(1, length(z1))
bg[grep('ETS', names(z1))] <- 'green'; cex[grep('ETS', names(z1))] <- 1
bg[grep('Stat', names(z1))] <- 'yellow'; cex[grep('Stat', names(z1))] <- 1
par(mar = c(10, 20, 10, 10)); plot(z1, z2, pch = 21, bg = bg, col = 'black', cex = cex, xaxt = 'n', xlim = c(min(z1), max(z1) * 1.5), ylim = c(min(z2), max(z2) * 1.5), yaxt = 'n', xlab = '', ylab = '')
abline(v = 0, h = 0, col = 'blue', lty = 2, lwd = 2)

i <- bg == 'yellow'
library(wordcloud); textplot(z1[i], z2[i], gsub('(.+?)\\/.+', '\\1', names(z1)[i]), xpd = TRUE, new = FALSE, adj = 0, cex = 1.5, col = 'black')

par(mar = c(10, 20, 10, 10)); plot(z1, z2, pch = 21, bg = bg, col = 'black', cex = 1, xaxt = 'n', xlim = c(min(z1), 0), ylim = c(min(z2), 0), yaxt = 'n', xlab = '', ylab = '')
abline(v = 0, h = 0, col = 'blue', lty = 2, lwd = 2)
i <- z1 < -0.01 & z2 < -0.01 & bg == 'gray'
library(wordcloud); textplot(z1[i], z2[i], gsub('(.+?)\\/.+', '\\1', names(z1)[i]), xpd = TRUE, new = FALSE, adj = 1, cex = 1.25, col = 'red')


# [2017-11-27] Generate the mononucleosome and NFR signals from the ATAC-seq data
project <- 'Etv2_MEFs'
project.dir <- sprintf('%s/projects/%s', Sys.getenv('SHARED'), project)
bam.files <- c(
#	'D2_EB_DOX' = sprintf('%s/mESCs_bams/d3_EB_Dox_1_S3.rmdup.bam', project.dir),
#	'D3_EB_DOX' = sprintf('%s/mESCs_bams/D3EB_F_P_6h_Dox-1_S12.rmdup.bam',project.dir),
#	'D3_EB' 		= sprintf('%s/mESCs_bams/D3EB_F_P_no_Dox-1_S10.rmdup.bam',project.dir),
#	'D2_EB' 		= sprintf('%s/mESCs_bams/d3_EB_no_dox_1_S1.rmdup.bam', project.dir),
#	'ESC_DOX' 	= sprintf('%s/mESCs_bams/ES_6h_Dox-1_S8.rmdup.bam', project.dir),
#	'ESC' 			= sprintf('%s/mESCs_bams/ES_no_Dox-1_S6.rmdup.bam', project.dir),
#	'D2_EB_DOX' = sprintf('%s/mESCs_bams/d3_EB_Dox_2_S4.rmdup.bam', project.dir),
#	'D3_EB_DOX' = sprintf('%s/mESCs_bams/D3EB_F_P_6h_Dox-2_S13.rmdup.bam', project.dir),
#	'D3_EB' 		= sprintf('%s/mESCs_bams/D3EB_F_P_no_Dox-2_S11.rmdup.bam', project.dir),
#	'D2_EB' 		= sprintf('%s/mESCs_bams/d3_EB_no_dox_2_S2.rmdup.bam', project.dir),
#	'ESC_DOX' 	= sprintf('%s/mESCs_bams/ES_6h_Dox-2_S9.rmdup.bam', project.dir),
#	'ESC' 			= sprintf('%s/mESCs_bams/ES_no_Dox-2_S7.rmdup.bam', project.dir),
	'MEF_DOX' 	= sprintf('%s/lab_and_public_MEFs_bams/3hDox1.rmdup.bam', project.dir),
	'MEF_DOX' 	= sprintf('%s/lab_and_public_MEFs_bams/3hDox2.rmdup.bam', project.dir),
	'MEF' 			= sprintf('%s/lab_and_public_MEFs_bams/noDox1.rmdup.bam', project.dir),
	'MEF' 			= sprintf('%s/lab_and_public_MEFs_bams/noDox2.rmdup.bam', project.dir)
#	'MEF_SRR1930168' = sprintf('%s/lab_and_public_MEFs_bams/SRR1930168.rmdup.bam', project.dir)
#	'Tie2pos' = '/panfs/roc/scratch/gongx030/ncbi/sra/ERR10/ERR1079/ERR1079314/ERR1079314.dedup.bam' # Tie2+ from E12.5
)
d <- data.frame(bam.file = bam.files, group = names(bam.files))
d <- transform(d, 
	coverage.file = gsub('.bam$', '.bed', bam.file), 
	coverage.mononuc.file = gsub('.bam$', 'mononuc.bed', bam.file),
	bam.index.file = sprintf('%s.bai', bam.file),
	bw.mononuc.file = gsub('.bam$', '.mononuc.bw', bam.file)
)
table(file.exists(d[, 'bam.file'])) # make sure all files exist

# [2017-11-27] find a subset of mononucleosome that are disappear in Dox induced samples
mef.mononuc.peaks.file <- sprintf('%s/MEF.mononuc.bed', project.dir)	# the mononuc positions in WT MEF
mef.mononuc.file <- sprintf('%s/MEF.mononuc.bw', project.dir)	# bigwig signal of mononucleosomes in MEF
mef_dox.mononuc.file <- sprintf('%s/MEF_DOX.mononuc.bw', project.dir)	# bigwig signal of mononucleosomes in MEF_DOX
mononuc_reduced_in_MEF_DOX.file <- sprintf('%s/mononuc_reduced_in_MEF_DOX.bed', project.dir)	# mononucleosomes reduced in Etv2 induced MEF compared with WT MEF

# [2017-11-27] mononucleosome peaks that are present in WT MEF
#source('chipseq.r'); cvg <- atac.coverage(d[, 'bam.file'], type = 'mononuc', genome = 'mm10')	# coverage of mononucleosome in WT MEF
source('chipseq.r'); cvg <- atac.coverage(d[d[, 'group'] == 'MEF', 'bam.file'], type = 'nfr', genome = 'mm10')	# coverage of NFR in WT MEF
#peaks <- slice(cvg, lower = 10)	# the mononucleosome positions at MEF
peaks <- slice(cvg, lower = 10)	# the mononucleosome positions at MEF
peaks <- reduce(as(peaks, 'GRanges'))
peaks <- add.seqinfo(peaks, genome = 'mm10')
peaks <- trim(resize(peaks, width = 250, fix = 'center'))# 250 bp regions surrounding the mononucleosomes
write.table(as.data.frame(peaks)[, 1:3], mef.mononuc.peaks.file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

# [2017-11-27] read the read count matrix for mononucleosome intervals in WT MEF (mef.mononuc.peaks.file)
peaks <- read.table(mef.mononuc.peaks.file, header = FALSE, sep = '\t')
peaks <- GRanges(seqnames = peaks[, 1], range = IRanges(peaks[, 2], peaks[, 3]))
source('aux.r'); peaks <- add.seqinfo(peaks, genome = 'mm10')
source('chipseq.r')
library(parallel); X <- do.call('cbind', mclapply(1:nrow(d), function(i){
	cvg <- atac.coverage(d[i, 'bam.file'], type = 'mononuc', genome = 'mm10')  # coverage of mononucleosome in WT MEF
	sum(cvg[peaks])
}, mc.cores = 2))
rownames(X) <- sprintf('%s:%d-%d', seqnames(peaks), start(peaks), end(peaks))

# [2017-11-27] Find the mononucleosome peaks that are significantly reduced in the MEF_DOX samples
library(DESeq2); dds <- DESeqDataSetFromMatrix(countData = X, colData = data.frame(group = d[, 'group']), design = ~ group)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <-	as.data.frame(results(dds, contrast = c('group', 'MEF_DOX', 'MEF')))
#rmn <- !is.na(res[, 'pvalue']) & res[, 'pvalue'] < 0.05 & res[, 'log2FoldChange'] < log2(1/1.5)
rmn <- !is.na(res[, 'pvalue']) & res[, 'pvalue'] < 0.05 & res[, 'log2FoldChange'] > log2(1.5)
write.table(as.data.frame(peaks)[rmn, 1:3], mononuc_reduced_in_MEF_DOX.file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)


# [2017-11-27] Find the positions of TFBS on reduced mononucleosome positions
library(chromVARmotifs) # https://github.com/GreenleafLab/chromVARmotifs
library(chromVAR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)
library(BiocParallel)
library(SummarizedExperiment)
source('filenames.r')
source('chipseq.r')
register(MulticoreParam(4)) # Use 8 cores
motif.set <- 'homer_pwms'
data(list = motif.set, package = 'chromVARmotifs')

# [2017-11-27] read the peaks, as well as the read count matrix associated with reduced mononuc positions
peaks <- read.table(mononuc_reduced_in_MEF_DOX.file, header = FALSE, sep = '\t')
peaks <- GRanges(seqnames = peaks[, 1], range = IRanges(peaks[, 2], peaks[, 3]))
motif_ix <- matchMotifs(get(motif.set), peaks, genome = 'mm10', out = 'positions')
library(parallel); X <- do.call('cbind', mclapply(1:nrow(d), function(i){
	cvg <- atac.coverage(d[i, 'bam.file'], type = 'mononuc', genome = 'mm10')  # coverage of mononucleosome in WT MEF
	sum(cvg[peaks])
}, mc.cores = 2))
rownames(X) <- sprintf('%s:%d-%d', seqnames(peaks), start(peaks), end(peaks))

# [2017-11-27] Find TFs associated with dynamic mononuc regions
rse <- SummarizedExperiment(assays = SimpleList(counts = X), rowRanges = peaks, colData = data.frame(group = d[, 'group']))
rse <- addGCBias(rse, genome = BSgenome.Mmusculus.UCSC.mm10)
dev <- computeDeviations(object = rse, annotations = matchMotifs(get(motif.set), peaks, genome = 'mm10'))
v <- computeVariability(dev)
#source('aux.r'); peaks <- add.seqinfo(peaks, genome = 'mm10')	# will cause error when calling matchMotifs
#js <- unique(c(order(v[, 'p_value'])[1:70], 82))
js <- unique(c(which(v[, 'p_value'] < 0.01), 82))

# the TF binding pattern
Z <- do.call('rbind', mclapply(js, function(j){
	if (length(motif_ix[[j]]) > 0){
		bed.file <- tempfile(); bed2.file <- tempfile(); bw.file <- tempfile(); bed3.file <- tempfile()
		write.table(cbind(as.data.frame(reduce(motif_ix[[j]], ignore.strand = TRUE))[, 1:3], 1), bed.file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
		command <- sprintf('sort -k1,1 -k2,2n %s > %s', bed.file, bed2.file)
		system(command)
		command <- sprintf('bedGraphToBigWig %s %s %s', bed2.file, genome.file('mm10', extend = TRUE), bw.file)
		system(command)
		command <- sprintf('bedtools intersect -a %s -b %s -wa > %s', mononuc_reduced_in_MEF_DOX.file, bed2.file, bed3.file)
		system(command)
		X <- bwtool.matrix(bed3.file, bw.file, up = 125, down = 125)
		X[is.na(X)] <- 0
		colSums(X)
	}else
		rep(0, 250)
}, mc.cores = 8))
Z[is.na(Z)] <- 0
rownames(Z) <- v[js, 'name']
Z <- t(apply(Z, 1, function(z) (z -min(z)) / (max(z) - min(z))))

library(pheatmap); library(gplots)
col <- colorpanel(100, low = 'blue', mid = 'black', high = 'yellow')
breaks <- seq(0, 1, length.out = 101)
pheatmap(Z, cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = FALSE, breaks = breaks, color = col, cellheight = 12, cellwidth = 1, fontsize_col = 15, fontsize_row = 12, show_colnames = FALSE, border_color = 'black', annotation_names_row = FALSE, annotation_names_col = FALSE)


# --------------------------------------------------------------------------
# [2017-11-28] Look at the distribution of Etv2/CTCF binding sites in dynamic changed NFR regions
# --------------------------------------------------------------------------
project <- 'Etv2_MEFs'
project.dir <- sprintf('%s/projects/%s', Sys.getenv('SHARED'), project)
bam.files <- c(
	'D2_EB_DOX' = sprintf('%s/mESCs_bams/d3_EB_Dox_1_S3.rmdup.bam', project.dir),
	'D3_EB_DOX' = sprintf('%s/mESCs_bams/D3EB_F_P_6h_Dox-1_S12.rmdup.bam',project.dir),
	'D3_EB' 		= sprintf('%s/mESCs_bams/D3EB_F_P_no_Dox-1_S10.rmdup.bam',project.dir),
	'D2_EB' 		= sprintf('%s/mESCs_bams/d3_EB_no_dox_1_S1.rmdup.bam', project.dir),
	'ESC_DOX' 	= sprintf('%s/mESCs_bams/ES_6h_Dox-1_S8.rmdup.bam', project.dir),
	'ESC' 			= sprintf('%s/mESCs_bams/ES_no_Dox-1_S6.rmdup.bam', project.dir),
	'D2_EB_DOX' = sprintf('%s/mESCs_bams/d3_EB_Dox_2_S4.rmdup.bam', project.dir),
	'D3_EB_DOX' = sprintf('%s/mESCs_bams/D3EB_F_P_6h_Dox-2_S13.rmdup.bam', project.dir),
	'D3_EB' 		= sprintf('%s/mESCs_bams/D3EB_F_P_no_Dox-2_S11.rmdup.bam', project.dir),
	'D2_EB' 		= sprintf('%s/mESCs_bams/d3_EB_no_dox_2_S2.rmdup.bam', project.dir),
	'ESC_DOX' 	= sprintf('%s/mESCs_bams/ES_6h_Dox-2_S9.rmdup.bam', project.dir),
	'ESC' 			= sprintf('%s/mESCs_bams/ES_no_Dox-2_S7.rmdup.bam', project.dir),
	'MEF_DOX' 	= sprintf('%s/lab_and_public_MEFs_bams/3hDox1.rmdup.bam', project.dir),
	'MEF_DOX' 	= sprintf('%s/lab_and_public_MEFs_bams/3hDox2.rmdup.bam', project.dir),
	'MEF' 			= sprintf('%s/lab_and_public_MEFs_bams/noDox1.rmdup.bam', project.dir),
	'MEF' 			= sprintf('%s/lab_and_public_MEFs_bams/noDox2.rmdup.bam', project.dir)
)
d <- data.frame(bam.file = bam.files, group = names(bam.files))
table(file.exists(d[, 'bam.file'])) # make sure all files exist

source('chipseq.r'); cvg <- mclapply(1:nrow(d), function(i) atac.coverage(d[i, 'bam.file'], type = 'nfr', genome = 'mm10'), mc.cores = 4)	# extract the NFR reads
#source('chipseq.r'); cvg <- mclapply(1:nrow(d), function(i) atac.coverage(d[i, 'bam.file'], type = 'mononuc', genome = 'mm10'), mc.cores = 4)	# extract the NFR reads
peaks <- slice(Reduce('+', cvg), lower = 10) # regions where merged reads > 10
peaks <- reduce(as(peaks, 'GRanges'))
peaks <- add.seqinfo(peaks, genome = 'mm10')
peaks <- trim(resize(peaks, width = 500, fix = 'center'))# 500 bp regions surrounding the mononucleosomes

# [2017-11-28] find the dynamic NFR using DESeq
X <- do.call('cbind', lapply(cvg, function(x) sum(x[peaks])))	# converting the Coverage object to a the read count matrix
rownames(X) <- sprintf('%s:%d-%d', seqnames(peaks), start(peaks), end(peaks))
library(DESeq2); dds <- DESeqDataSetFromMatrix(countData = X, colData = data.frame(group = d[, 'group']), design = ~ group)
dds <- DESeq(dds)
res <-	as.data.frame(results(dds, contrast = c('group', 'MEF_DOX', 'MEF')))	
dr <- !is.na(res[, 'pvalue']) & res[, 'pvalue'] < 0.001 & abs(res[, 'log2FoldChange']) > log2(1.5)	# dynamically changed NFR


# [2017-11-28] the density matrix for the dynamic NFR, using a wider range (2000 bp) for better visualization
peaks2 <- trim(resize(peaks, width = 2000, fix = 'center'))
Y <- t(sapply(Reduce('+', cvg[d[, 'group'] == 'MEF'])[peaks2[dr]], as.vector))
Y.dox <- t(sapply(Reduce('+', cvg[d[, 'group'] == 'MEF_DOX'])[peaks2[dr]], as.vector))
hh <- order(res[dr, 'log2FoldChange'])	# sort the peaks based on the fold change at the center 500bp region
dY <- Y.dox - Y
dY <- t(scale(t(dY)))
breaks <- c(-10, seq(-5, 5, length.out = 99), 10)
image(t(dY[hh, ]), col = colorpanel(100, low = 'blue', mid = 'black', high = 'yellow'), axes = FALSE)


# [2017-11-27] Find the positions of TFBS on reduced mononucleosome positions
library(chromVARmotifs) # https://github.com/GreenleafLab/chromVARmotifs
library(chromVAR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)
library(BiocParallel)
library(SummarizedExperiment)
source('filenames.r')
source('chipseq.r')
register(MulticoreParam(4)) # Use 8 cores
motif.set <- 'homer_pwms'
data(list = motif.set, package = 'chromVARmotifs')

# [2017-11-27] First Find TFs associated with dynamic intervals
rse <- SummarizedExperiment(assays = SimpleList(counts = X[dr, ]), rowRanges = peaks[dr], colData = data.frame(group = d[, 'group']))
rse <- addGCBias(rse, genome = BSgenome.Mmusculus.UCSC.mm10)
dev <- computeDeviations(object = rse, annotations = matchMotifs(get(motif.set), peaks[dr], genome = 'mm10'))
v <- computeVariability(dev)

# !!!!!!!!!!!!!!!!!!!!!!!!!!!
#save.image(file = sprintf('%s/R_sessions/%s.RData', project.dir, session.name))
#load(sprintf('%s/R_sessions/%s.RData', project.dir, session.name))
# !!!!!!!!!!!!!!!!!!!!!!!!!!!

h.high.in.dox <- !is.na(res[, 'pvalue']) & res[, 'pvalue'] < 0.001 & res[, 'log2FoldChange'] > log2(1.5)	# NFR high in DOX
h.high.in.wt <- !is.na(res[, 'pvalue']) & res[, 'pvalue'] < 0.001 & res[, 'log2FoldChange'] < log2(1 / 1.5)	# NFR high in WT

#j <- 'Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer(0.967)'
#j <- 'FOXA1:AR(Forkhead,NR)/LNCAP-AR-ChIP-Seq(GSE27824)/Homer'
j <- 'Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer'
#j <- 'CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer'
#j <- 'ETS(ETS)/Promoter/Homer'
#j <- 'NeuroD1(bHLH)/Islet-NeuroD1-ChIP-Seq(GSE30298)/Homer'
#j <- 'GATA(Zf),IR4/iTreg-Gata3-ChIP-Seq(GSE20898)/Homer'
#j <- 'STAT5(Stat)/mCD4+-Stat5-ChIP-Seq(GSE12346)/Homer'
#j <- 'TATA-Box(TBP)/Promoter/Homer'
#j <- 'PAX3:FKHR-fusion(Paired,Homeobox)/Rh4-PAX3:FKHR-ChIP-Seq(GSE19063)/Homer'
#j <- 'Zic(Zf)/Cerebellum-ZIC1.2-ChIP-Seq(GSE60731)/Homer'

Z2 <- do.call('rbind', lapply(list(h.high.in.dox, h.high.in.wt), function(h){
	mm <- matchMotifs(get(motif.set)[j], peaks2[h], genome = 'mm10', out = 'positions')[[1]]
	cvg.motif <- coverage(add.seqinfo(mm, genome = 'mm10'), weight = values(mm)[['score']]) # it is import to add the seqinfo; by default, motif_ix returned from matchMotifs has no seqinfo
	Z <- as(as(cvg.motif[peaks2[h]], 'RleViews'), 'matrix')
	Z <- t(apply(Z, 1, function(z) (z - min(z)) / (max(z) - min(z))))
	Z[is.na(Z)] <- 0
	colMeans(Z)
}))

plot(1:ncol(Z2), ylim = c(0, max(Z2)), xaxt = 'n', yaxt = 'n')
lines(Z2[1, ], lwd = 4, col = 'blue')
lines(Z2[2, ], lwd = 4, col = 'red')
abline(v = c(750, 1250), lwd = 3, lty = 2)


# [2017-11-30] Motifmatch results from mouse DHS 300 sites, with the positional information
project <- 'Etv2_MEFs'
project.dir <- sprintf('%s/projects/%s', Sys.getenv('SHARED'), project)
mm.file <- sprintf('%s/matchMotifs_mouse_dhs.300.RData', project.dir); genome <- 'mm10'; motif.set <- 'homer_pwms'
dhs.bed.file <- sprintf('%s/projects/sc_atac/mouse_dhs.300.bed', Sys.getenv('SHARED'))
source('chipseq.r'); peaks <- read.table(dhs.bed.file, header = FALSE, sep = '\t')
peaks <- GRanges(seqnames = peaks[, 1], range = IRanges(peaks[, 2], peaks[, 3]))
peaks <- add.seqinfo(peaks, genome = genome)
mm <- matchMotifs(get(motif.set), peaks, genome = 'mm10', out = 'positions')
save(mm, file = mm.file)


# [2017-12-01] MNase coverage at DHS 300 sites in MEF and ESC
project <- 'Etv2_MEFs'; session.name <- 'Dox_specific_TF_at_nucleosome'; 
project.dir <- sprintf('%s/projects/%s', Sys.getenv('SHARED'), project)
genome <- 'mm10'; mnase.width <- 500
dhs.bed.file <- sprintf('%s/projects/sc_atac/mouse_dhs.300.bed', Sys.getenv('SHARED'))
mnase.coverage.files <- c(
	'MEF' = sprintf('%s/projects/MNase/MEF_Coverage_pileup.bw', Sys.getenv('SHARED')),
	'ESC' = sprintf('%s/projects/MNase/mESCs_Coverage_pileup.bw', Sys.getenv('SHARED'))
)
source('chipseq.r'); peaks <- read.table(dhs.bed.file, header = FALSE, sep = '\t')
peaks <- GRanges(seqnames = peaks[, 1], range = IRanges(peaks[, 2], peaks[, 3]))
peaks <- add.seqinfo(peaks, genome = genome)
cell <- 'MEF'
source('chipseq.r'); cvg.mnase <- import(mnase.coverage.files[cell], which = trim(reduce(resize(peaks, width = mnase.width , fix = 'center'))), as = 'RleList')  # the depth coverage of MNase-seq
mnase.cvg.file <- sprintf('%s/mnase_coverage_%s.RData', project.dir, cell)
save(cvg.mnase, file = mnase.cvg.file)


# [2017-11-29] For ATAC-seq peaks, identify the distance of each TF motif to nucleosome positions
project <- 'Etv2_MEFs'; session.name <- 'Dox_specific_TF_at_nucleosome'; 
project.dir <- sprintf('%s/projects/%s', Sys.getenv('SHARED'), project)
genome <- 'mm10'; atac.width <- 250; min.depth <- 20; mnase.width <- 100; motif.set <- 'homer_pwms'
bam.files <- c(
	'D2_EB_DOX' = sprintf('%s/mESCs_bams/d3_EB_Dox_1_S3.rmdup.bam', project.dir),
	'D3_EB_DOX' = sprintf('%s/mESCs_bams/D3EB_F_P_6h_Dox-1_S12.rmdup.bam',project.dir),
	'D3_EB' 		= sprintf('%s/mESCs_bams/D3EB_F_P_no_Dox-1_S10.rmdup.bam',project.dir),
	'D2_EB' 		= sprintf('%s/mESCs_bams/d3_EB_no_dox_1_S1.rmdup.bam', project.dir),
	'ESC_DOX' 	= sprintf('%s/mESCs_bams/ES_6h_Dox-1_S8.rmdup.bam', project.dir),
	'ESC' 			= sprintf('%s/mESCs_bams/ES_no_Dox-1_S6.rmdup.bam', project.dir),
	'D2_EB_DOX' = sprintf('%s/mESCs_bams/d3_EB_Dox_2_S4.rmdup.bam', project.dir),
	'D3_EB_DOX' = sprintf('%s/mESCs_bams/D3EB_F_P_6h_Dox-2_S13.rmdup.bam', project.dir),
	'D3_EB' 		= sprintf('%s/mESCs_bams/D3EB_F_P_no_Dox-2_S11.rmdup.bam', project.dir),
	'D2_EB' 		= sprintf('%s/mESCs_bams/d3_EB_no_dox_2_S2.rmdup.bam', project.dir),
	'ESC_DOX' 	= sprintf('%s/mESCs_bams/ES_6h_Dox-2_S9.rmdup.bam', project.dir),
	'ESC' 			= sprintf('%s/mESCs_bams/ES_no_Dox-2_S7.rmdup.bam', project.dir),
	'MEF_DOX' 	= sprintf('%s/lab_and_public_MEFs_bams/3hDox1.rmdup.bam', project.dir),
	'MEF_DOX' 	= sprintf('%s/lab_and_public_MEFs_bams/3hDox2.rmdup.bam', project.dir),
	'MEF' 			= sprintf('%s/lab_and_public_MEFs_bams/noDox1.rmdup.bam', project.dir),
	'MEF' 			= sprintf('%s/lab_and_public_MEFs_bams/noDox2.rmdup.bam', project.dir)
)
d <- data.frame(group = names(bam.files), bam.file = bam.files)
group2cell <- c('MEF_DOX' = 'MEF', 'MEF' = 'MEF', 'D2_EB_DOX' = 'ESC', 'D2_EB' = 'ESC', 'ESC_DOX' = 'ESC', 'ESC' = 'ESC', 'D3_EB' = 'ESC', 'D3_EB_DOX' = 'ESC')
mnase.coverage.files <- c(
	'MEF' = sprintf('%s/projects/MNase/MEF_Coverage_pileup.bw', Sys.getenv('SHARED')),
	'ESC' = sprintf('%s/projects/MNase/mESCs_Coverage_pileup.bw', Sys.getenv('SHARED'))
)
d <- transform(d, cell = group2cell[group])
d <- transform(d, mnase.coverage.file = mnase.coverage.files[cell])
d <- transform(d, nucdist.file = gsub('bam$', 'nucdist.RData', bam.file))
table(file.exists(d[, 'bam.file'])) # make sure all files exist
peak.file <- sprintf('%s/ATAC_peaks_min.depth=%d_atac.width=%d_version=20171205a.RData', project.dir, min.depth, atac.width)	# the merged ATAC-seq peak union
matchMotifs.file <- gsub('ATAC_peaks', 'matchMotifs', peak.file)	# the matchMotifs results for all the union ATAC-seq peak set


# [2017-12-05] Determine the ATAC peaks and save the read count matrix
source('chipseq.r'); se <- atac.read.peaks(d[, 'bam.file'], genome = genome, min.depth = min.depth, atac.width = atac.width)
se <- add.seqinfo(se, genome)
save(se, file = peak.file)


# [2017-12-04] Pre-Compute the moatchMotifs results (with positions) on ATAC-seq peaks
se <- local({get(load(peak.file))})
mm <- matchMotifs(get(motif.set), se, genome = genome, out = 'positions')
save(mm, file = matchMotifs.file)

# [2017-12-05] Get the Mnase coverage for all ATAC-seq peaks at ESC and MEF
mm <- local({get(load(matchMotifs.file))})
cell <- 'MEF'; treatment <- 'MEF_DOX'; control <- 'MEF'; min.depth <- 10; atac.width <- 50; mnase.width <- 500
#cell <- 'ESC'; treatment <- 'D2_EB_DOX'; control <- 'D2_EB'; min.depth <- 10; atac.width <- 50; mnase.width <- 500
m <- d[, 'group'] %in% c(treatment, control)
source('chipseq.r'); cvg <- mclapply(which(m), function(i) atac.coverage(d[i, 'bam.file'], genome = 'mm10'), mc.cores = 4)	# extract the NFR reads

library(DESeq2)
mclapply(1:length(mm), function(h){
	r <- add.seqinfo(mm[[h]], genome = 'mm10')
	X <- do.call('cbind', lapply(cvg, function(c) sum(c[r])))	# the read count matrix on TFBS
	dds <- DESeqDataSetFromMatrix(countData = X, colData = d[m, 'group', drop = FALSE], design = ~ group)
	dds <- estimateSizeFactors(dds)
	dds <- DESeq(dds)
	res <-	as.data.frame(results(dds, contrast = c('group', treatment, control)))
	r <- trim(resize(r, width = mnase.width, fix = 'center'))
	cvg.mnase <- import(mnase.coverage.files[cell], which = r, as = 'RleList')  # the depth coverage of MNase-seq
	Y <- as(as(cvg.mnase[r], 'RleViews'), 'matrix')
	values(r)[['pvalue']] <- res[, 'pvalue']
	values(r)[['log2FoldChange']] <- res[, 'log2FoldChange']
	values(r)[['mnase']] <- Y
	se <- SummarizedExperiment(assays = list(mnase = X), rowRanges = r, colData = d[m, ])
	mnase.file <- sprintf('%s/atac_cell=%s_traetment=%s_control=%s_tf=%d.RData', project.dir, cell, treatment, control, h)
	cat(sprintf('[%s] writing %s\n', Sys.time(), mnase.file))
	save(se, file = mnase.file)
}, mc.cores = 2)


# [2017-12-06] Demonstrating that Etv2 peaks with 

mm <- local({get(load(matchMotifs.file))})
cell <- 'MEF'; treatment <- 'MEF_DOX'; control <- 'MEF'
#cell <- 'ESC'; treatment <- 'D2_EB_DOX'; control <- 'D2_EB';
tf <- 'Oct4(POU,Homeobox)/mES-Oct4-ChIP-Seq(GSE11431)/Homer'
h <- which(names(mm) == tf)
mnase.file <- sprintf('%s/atac_cell=%s_traetment=%s_control=%s_tf=%d.RData', project.dir, cell, treatment, control, h)
source('chipseq.r'); se <- local({get(load(mnase.file))})
Y <- values(rowRanges(se))[['mnase']]
Y <- t(scale(t(log(Y + 1)))); Y[is.na(Y)] <- 0
Y <- (Y + Y[, ncol(Y):1])[, round(ncol(Y) / 2):ncol(Y)]
up <- !is.na(values(rowRanges(se))[['pvalue']]) & values(rowRanges(se))[['pvalue']] < 0.05 & values(rowRanges(se))[['log2FoldChange']] > 0 # high in treatment
down <- !is.na(values(rowRanges(se))[['pvalue']]) & values(rowRanges(se))[['pvalue']] < 0.05 & values(rowRanges(se))[['log2FoldChange']] < 0 # high in control
xx <- rbind(up = colMeans(Y[up, ]), average = colMeans(Y), down = colMeans(Y[down, ]))
plot(NA, xlim = c(1, ncol(xx)), ylim = c(range(xx)))
cols <- c('up' = 'red', 'down' = 'blue', 'average' = 'black')
sapply(1:nrow(xx), function(i) lines(x = 1:ncol(xx), y = xx[i, ], col = cols[i]))


# [2017-12-06] Heatmap of TF ~ Mnase pattern on ATAC-seq peaks are significantly open on Etv2 induction
mm <- local({get(load(matchMotifs.file))})
cell <- 'MEF'; treatment <- 'MEF_DOX'; control <- 'MEF'
#cell <- 'ESC'; treatment <- 'D2_EB_DOX'; control <- 'D2_EB';
source('chipseq.r'); nd <- mclapply(1:length(mm), function(h){
	tryCatch({
		mnase.file <- sprintf('%s/atac_cell=%s_traetment=%s_control=%s_tf=%d.RData', project.dir, cell, treatment, control, h)
		cat(sprintf('[%s] reading %s\n', Sys.time(), mnase.file))
		se <- local({get(load(mnase.file))})
		up <- !is.na(values(rowRanges(se))[['pvalue']]) & values(rowRanges(se))[['pvalue']] < 0.05 & values(rowRanges(se))[['log2FoldChange']] > 0 # high in treatment
		down <- !is.na(values(rowRanges(se))[['pvalue']]) & values(rowRanges(se))[['pvalue']] < 0.05 & values(rowRanges(se))[['log2FoldChange']] < 0 # high in treatment
		Y <- values(rowRanges(se))[['mnase']]
		Y <- t(scale(t(log(Y + 1)))); Y[is.na(Y)] <- 0
		Y <- (Y + Y[, ncol(Y):1])[, round(ncol(Y) / 2):ncol(Y)]
		list(up = colMeans(Y[up, ]), down = colMeans(Y[down, ]), mean = colMeans(Y), n.up = sum(up), n.down = sum(down), n = nrow(Y))
	}, error = function(e) NULL)
}, mc.cores = 8, mc.preschedule = FALSE)
nd.file <- sprintf('%s/nucdist_cell=%s.RData', project.dir, cell)
names(nd) <- gsub('(.+?)/.+', '\\1', names(mm))
save(nd, file = nd.file)

# [2017-12-06] Look at the % of changed open TFBS
cell <- 'ESC'; treatment <- 'D2_EB_DOX'; control <- 'D2_EB';
nd.file <- sprintf('%s/nucdist_cell=%s.RData', project.dir, cell)
nd <- local({get(load(nd.file))})
Y <- sapply(nd, function(x) c(up = x$n.up / x$n, down = x$n.down / x$n))
i <- order(Y['up', ], decreasing = TRUE)[1:20]
dev.new(width = 7, height = 7); par(mar = c(5, 20, 2, 2)); barplot(rev(Y['up', i]), horiz = TRUE, las = 2, col = 'black', xpd = FALSE, xlim = c(0.03, 0.04))
i <- order(Y['down', ], decreasing = TRUE)[1:20]
dev.new(width = 7, height = 7); par(mar = c(5, 20, 2, 2)); barplot(rev(Y['down', i]), horiz = TRUE, las = 2, col = 'black', xpd = FALSE, xlim = c(0.06, 0.07))

Y2 <- Y[, !grepl('ETS', colnames(Y))]
i <- order(Y2['up', ], decreasing = TRUE)[1:20]
dev.new(width = 7, height = 7); par(mar = c(5, 20, 2, 2)); barplot(rev(Y2['up', i]), horiz = TRUE, las = 2, col = 'black', xpd = FALSE, xlim = c(0.03, 0.035))

cell <- 'MEF'; treatment <- 'MEF_DOX'; control <- 'MEF'
nd.file <- sprintf('%s/nucdist_cell=%s.RData', project.dir, cell)
nd <- local({get(load(nd.file))})
Y <- sapply(nd, function(x) c(up = x$n.up / x$n, down = x$n.down / x$n))
i <- order(Y['up', ], decreasing = TRUE)[1:20]
dev.new(width = 7, height = 7); par(mar = c(5, 20, 2, 2)); barplot(rev(Y['up', i]), horiz = TRUE, las = 2, col = 'black', xpd = FALSE, xlim = c(0.04, 0.065))
i <- order(Y['down', ], decreasing = TRUE)[1:20]
dev.new(width = 8, height = 7); par(mar = c(5, 25, 2, 2)); barplot(rev(Y['down', i]), horiz = TRUE, las = 2, col = 'black', xpd = FALSE, xlim = c(0.04, 0.06))

Y2 <- Y[, !grepl('ETS', colnames(Y))]
i <- order(Y2['up', ], decreasing = TRUE)[1:20]
dev.new(width = 7, height = 7); par(mar = c(5, 20, 2, 2)); barplot(rev(Y2['up', i]), horiz = TRUE, las = 2, col = 'black', xpd = FALSE, xlim = c(0.04, 0.045))

# most increased non-ESC factors


# [2017-12-06] Look at the commonly gain non-ETS TF sites
cell <- 'ESC'; treatment <- 'D2_EB_DOX'; control <- 'D2_EB';
nd.file <- sprintf('%s/nucdist_cell=%s.RData', project.dir, cell)
nd.esc <- local({get(load(nd.file))})
cell <- 'MEF'; treatment <- 'MEF_DOX'; control <- 'MEF'
nd.file <- sprintf('%s/nucdist_cell=%s.RData', project.dir, cell)
nd.mef <- local({get(load(nd.file))})
Y.esc <- sapply(nd.esc, function(x) c(up = x$n.up / x$n, down = x$n.down / x$n))
Y.mef <- sapply(nd.mef, function(x) c(up = x$n.up / x$n, down = x$n.down / x$n))

i <- which(!grepl('ETS', names(nd.mef)))
Y.esc <- Y.esc[, i]; Y.mef <- Y.mef[, i]

bg <- rep('black', ncol(Y.esc))
bg[grep('GATA', names(nd.mef))] <- 'red'
bg[grep('Forkhead', names(nd.mef))] <- 'yellow'
i <- bg == 'black'
plot(Y.esc['up', i], Y.mef['up', i], bg = bg[i], pch = 21, col = 'darkgray', xlim = range(Y.esc['up', ]), ylim = range(Y.mef['up', ]))
points(Y.esc['up', !i], Y.mef['up', !i], bg = bg[!i], pch = 21, col = 'darkgray')


#cell <- 'MEF'; treatment <- 'MEF_DOX'; control <- 'MEF'
cell <- 'ESC'; treatment <- 'D2_EB_DOX'; control <- 'D2_EB';
nd.file <- sprintf('%s/nucdist_cell=%s.RData', project.dir, cell)
nd <- local({get(load(nd.file))})
Z <- do.call('rbind', lapply(nd, function(x) x$up))
rn <- rownames(Z)
rownames(Z) <- sprintf('C%d', 1:nrow(Z))
library(pheatmap); library(gplots)
col <- colorpanel(100, low = 'blue', mid = 'black', high = 'yellow')
df.row <- data.frame(up = sapply(nd, function(x) x$n.up / x$n), down = sapply(nd, function(x) x$n.down / x$n))
rownames(df.row) <- rownames(Z)
breaks <- c(-10, seq(-0.5, 0.5, length.out = 99), 10)

j1 <- order(df.row[, 'up'], decreasing = TRUE)[1:25]
j2 <- order(df.row[, 'down'], decreasing = TRUE)[1:25]
j <- unique(c(j1, j2, c(43)))
pheatmap(Z[j, ], cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = FALSE, breaks = breaks, color = col, cellheight = 11, cellwidth = 1, fontsize_col = 9, fontsize_row = 11, show_colnames = FALSE, border_color = 'black', annotation_names_row = TRUE, annotation_names_col = FALSE, annotation_row = df.row, labels_row = rn[j])


# [2017-12-06] Look at the Etv2 sites relased on Etv2 induction
#cell <- 'MEF'; treatment <- 'MEF_DOX'; control <- 'MEF'
cell <- 'ESC'; treatment <- 'D2_EB_DOX'; control <- 'D2_EB';
h <- 82
h <- 102 #GATA
mnase.file <- sprintf('%s/atac_cell=%s_traetment=%s_control=%s_tf=%d.RData', project.dir, cell, treatment, control, h)
source('chipseq.r'); se <- local({get(load(mnase.file))})
Y <- values(rowRanges(se))[['mnase']]
Y <- t(scale(t(Y))); Y[is.na(Y)] <- 0
Y <- (Y + Y[, ncol(Y):1])[, round(ncol(Y) / 2):ncol(Y)]
fc <- values(rowRanges(se))[['log2FoldChange']]; Y <- Y[!is.na(fc), ]; fc <- fc[!is.na(fc)]
sp <- split(1:nrow(Y), list(cut(fc, quantile(fc, seq(0, 1, by = 0.02)))))
Y2 <- do.call('rbind', lapply(sp, function(i) colMeans(Y[i, ])))
rownames(Y2) <- sprintf('C%d', 1:nrow(Y2))
df.row <- data.frame(fc = sapply(sp, function(i) mean(fc[i])))
rownames(df.row) <- rownames(Y2)
col <- colorpanel(100, low = 'blue', mid = 'black', high = 'yellow')
breaks <- c(-20, seq(-0.5, 0.5, length.out = 99), 20)
pheatmap(Y2[order(df.row[, 'fc']), ], breaks = breaks, cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = FALSE, color = col, cellheight = 8, cellwidth = 1, fontsize_col = 8, show_colnames = FALSE, border_color = 'black', annotation_names_row = FALSE, annotation_names_col = FALSE, annotation_row = df.row)

# [2017-12-06] Look for Etv2 cofactors that specifically open the chomosomes 
#cell <- 'MEF'; treatment <- 'MEF_DOX'; control <- 'MEF'
cell <- 'ESC'; treatment <- 'D2_EB_DOX'; control <- 'D2_EB';
h <- 82	# Etv2
mnase.file <- sprintf('%s/atac_cell=%s_traetment=%s_control=%s_tf=%d.RData', project.dir, cell, treatment, control, h)
source('chipseq.r'); se <- local({get(load(mnase.file))})
up <- !is.na(values(rowRanges(se))[['pvalue']]) & values(rowRanges(se))[['pvalue']] < 0.05 & values(rowRanges(se))[['log2FoldChange']] > 0 # high in treatment
r <- resize(rowRanges(se), width = 50, fix = 'center')
mm <- matchMotifs(get(motif.set), r, genome = genome)
B <- assay(mm)
Y <- values(rowRanges(se))[['mnase']]
Y <- t(scale(t(Y))); Y[is.na(Y)] <- 0
Y <- (Y + Y[, ncol(Y):1])[, round(ncol(Y) / 2):ncol(Y)]
fc <- values(rowRanges(se))[['log2FoldChange']]
B <- B[!is.na(fc), ]; Y <- Y[!is.na(fc), ]; up <- up[!is.na(fc)]; fc <- fc[!is.na(fc)]; 

#res <- do.call('rbind', mclapply(1:ncol(B), function(i) data.frame(name = colnames(B)[i], n = sum(B[, i]), pvalue = wilcox.test(fc[B[, i] & up], fc[B[, i] & !up], alternative = 'greater')$p.value)))
res <- do.call('rbind', mclapply(1:ncol(B), function(i) data.frame(name = colnames(B)[i], n = sum(B[, i]), pvalue = fisher.test(table(factor(up, c(FALSE, TRUE)), factor(B[, i] > 0, c(FALSE, TRUE))), alternative = 'greater')$p.value), mc.cores = 8))
res <- res[order(res[, 'pvalue']), ]
res <- res[!grepl('ETS', res[, 'name']), ]
x <- -log(res[1:20, 'pvalue'])
names(x) <- gsub('(.+?)/.+', '\\1', res[1:20, 'name'])
dev.new(width = 7, height = 7); par(mar = c(5, 20, 2, 2)); barplot(rev(x), horiz = TRUE, las = 2, col = 'black', xpd = FALSE)



Z.up <- do.call('rbind', mclapply(1:ncol(B), function(i) colMeans(Y[B[, i] & up, , drop = FALSE]) - colMeans(Y[B[, i], ]), mc.cores = 8))
Z.mean <- do.call('rbind', mclapply(1:ncol(B), function(i) colMeans(Y[B[, i], ]), mc.cores = 8))

col <- colorpanel(100, low = 'blue', mid = 'black', high = 'yellow')
breaks <- c(-20, seq(-1, 1, length.out = 99), 20)
Z <- Z.up - Z.mean
rownames(Z) <- gsub('(.+?)/.+', '\\1', names(get(motif.set)))
#i <- c(order(rowSums(Z[, 1:20]), decreasing = FALSE)[1:50], order(rowSums(Z[, 1:20]), decreasing = TRUE)[1:50])
i <- order(rowSums(Z[, 1:20]), decreasing = TRUE)[1:50]
Z2 <- Z[i, ]
pheatmap(Z2, breaks = breaks, cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = FALSE, color = col, cellheight = 10, cellwidth = 1, fontsize_col = 10, show_colnames = FALSE, border_color = 'black', annotation_names_row = FALSE, annotation_names_col = FALSE)

h <- 136
plot(colMeans(Y[B[, h] & up, ]))


h2 <- 150
Y2 <- Y[B[, h2], ]
rownames(Y2) <- sprintf('C%d', 1:nrow(Y2))
df.row <- data.frame(fc = fc[B[, h2]])
rownames(df.row) <- rownames(Y2)
pheatmap(Y2[order(df.row[, 'fc']), ], breaks = breaks, cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = FALSE, color = col, cellheight = 1000 / nrow(Y2), cellwidth = 1, fontsize_col = 8, show_colnames = FALSE, border_color = 'black', annotation_names_row = FALSE, annotation_names_col = FALSE, annotation_row = df.row)





#z0 <- sapply(1:1000, function(seed) mean(Y[sample(1:nrow(Y), sum(i)), 1:10]))
#z1 <- mean(Y[i, 1:10])

plot(colMeans(Y[sample(1:nrow(Y), sum(i)), ]) - colMeans(Y))

(Y[i, 1:10]) - mean(Y[sample(1:nrow(Y), sum(i)), 1:10])




values(se)[['pvalue']] <- res[, 'pvalue']


# [2017-12-04] Determine the peaks in ESC and MEF; including both dynamic and non-dynamic peaks (as long as the p-value is not NA)
treatment <- 'D2_EB_DOX'; control <- 'D2_EB'
#treatment <- 'MEF_DOX'; control <- 'MEF'
se <- local({get(load(peak.file))})
m <- d[, 'group'] %in% c(treatment, control)
values(se)[['log2FoldChange']] <- res[, 'log2FoldChange']
peak2.file <- gsub('.RData', sprintf('_treatment=%s_control=%s.RData', treatment, control), peak.file)
save(se, file = peak2.file)





peak2.file <- gsub('.RData', sprintf('_treatment=%s_control=%s.RData', treatment, control), peak.file)
se <- local({get(load(peak2.file))})
cvg.mnase <- import(mnase.coverage.files[cell], which = rowRanges(se), as = 'RleList')  # the depth coverage of MNase-seq

#tf <- 'Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer'
#tf <- 'Oct4(POU,Homeobox)/mES-Oct4-ChIP-Seq(GSE11431)/Homer'
tf <- 'Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer(0.967)'
source('chipseq.r'); p <- add.seqinfo(mm[[h]], genome = genome)
r <- trim(resize(p, width = 50, fix = 'center'))
Y <- as(as(cvg.mnase[r], 'RleViews'), 'matrix')
Y <- t(apply(Y, 1, function(y) (y - min(y)) / (max(y) - min(y)))); Y[is.na(Y)] <- 0

B <- as(findOverlaps(peaks, r), 'matrix')
B <- sparseMatrix(i = B[, 'queryHits'], j = B[, 'subjectHits'], dims = c(length(peaks), length(r)))
X <- as.matrix(t(B) %*% assay(peaks))

Z <- as.matrix(t(X %*% Diagonal(x = 1 / colSums(X))) %*% Y)
n.perm <- 200; Zm <- Reduce('+', mclapply(1:n.perm, function(i){
	Xp <- do.call('cbind', lapply(1:ncol(X), function(m) rmultinom(1, size = sum(X[, m]), prob = rep(1 / nrow(X), nrow(X))) / sum(X[, m])))
	t(Xp) %*%  Y
}, mc.cores = 4)) / n.perm

dZ <- t(scale(t(Z - Zm)))
#dZ <- t(scale(t(Zm)))
#dZ <- Zm
dZ <- (dZ + dZ[, ncol(dZ):1])[, round(ncol(dZ) / 2):ncol(dZ)]
image(t(as.matrix(dZ)), col = colorpanel(100, low = 'blue', mid = 'black', high = 'yellow'))

plot((Z - Zp)[1, ])


E <- rowSums(X) %o% colSums(X) / sum(X)








#tf <- 'Oct4(POU,Homeobox)/mES-Oct4-ChIP-Seq(GSE11431)/Homer'
#tf <- 'Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer'
tf <- 'Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer(0.967)'
#tf <- 'FOXA1:AR(Forkhead,NR)/LNCAP-AR-ChIP-Seq(GSE27824)/Homer'
h <- which(names(get(motif.set)) == tf)
p <- add.seqinfo(mm[[h]], genome = genome)
X <- do.call('cbind', lapply(cvg, function(x) sum(x[p])))	
r <- trim(resize(p, width = mnase.width , fix = 'center'))
Y <- as(as(cvg.mnase[r], 'RleViews'), 'matrix')
Y <- (Y + Y[, ncol(Y):1])[, round(ncol(Y) / 2):ncol(Y)]
Y <- do.call('rbind', lapply(1:nrow(d2), function(m) Matrix::colMeans(Diagonal(x = X[, m])) %*% Y ))
Y <- t(apply(Y, 1, function(y) (y - min(y)) / (max(y) - min(y))))
plot(colMeans(Y[1:2, ]))
plot(colMeans(Y[3:4, ]))

	nucdist.file <- sprintf('%s/%s_TF=%d.RData', project.dir, cell, h)
	cat(sprintf('[%s] writing %s\n', Sys.time(), nucdist.file))
	save(Y, file = nucdist.file)
#}


#	[2017-12-01] 
source('chipseq.r');
motif.set <- 'homer_pwms'
get(motif.set)
#tf <- 'FOXA1:AR(Forkhead,NR)/LNCAP-AR-ChIP-Seq(GSE27824)/Homer'
#tf <- 'Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer'
#tf <- 'OCT4-SOX2-TCF-NANOG(POU,Homeobox,HMG)/mES-Oct4-ChIP-Seq(GSE11431)/Homer'
tf <- 'Oct4(POU,Homeobox)/mES-Oct4-ChIP-Seq(GSE11431)/Homer'
#tf <- 'CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer'
#tf <- 'Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer(0.967)'
#tf <- 'bHLHE40(bHLH)/HepG2-BHLHE40-ChIP-Seq(GSE31477)/Homer'
h <- which(names(get(motif.set)) == tf)

Y <- do.call('rbind', lapply(c('MEF', 'ESC'), function(cell){
	nucdist.file <- sprintf('%s/%s_TF=%d.RData', project.dir, cell, h)
	y <- local({get(load(nucdist.file))})
	rownames(y) <- d[d[, 'cell'] == cell, 'group']
	y
}))
Y <- t(apply(Y, 1, function(y) (y - min(y)) / (max(y) - min(y))))
Y <- aggregate(Y, list(rownames(Y)), mean)
rownames(Y) <- Y[, 1]; Y <- as.matrix(Y[, -1])

Z <- rbind(Y['MEF', ] - Y['MEF_DOX', ], Y['D2_EB', ] - Y['D2_EB_DOX', ]); rownames(Z) <- c('MEF', 'ESC')
dev.new(width = 7, height = 5); plot(NA, xlim = c(1, ncol(Y)), ylim = range(Z), xlab = '', ylab = '')
col <- c('MEF' = 'blue', 'ESC' = 'red')
lapply(1:nrow(Z), function(i) lines(1:ncol(Z), Z[i, ], col = col[rownames(Z)[i]], lty = 1, lwd = 4))
abline(v = c(75, 150), h = 0, lty = 2, col = 'black', lwd = 3)

plot(NA, xlim = c(1, ncol(Y)), ylim = range(Y[1:2, ]), xlab = '', ylab = '')
lines(1:ncol(Y), Y[1, ], col = 'black', lty = 1, lwd = 4)
lines(1:ncol(Y), Y[2, ], col = 'black', lty = 2, lwd = 4)

plot(NA, xlim = c(1, ncol(Z)), ylim = range(Z[2, ]), xlab = '', ylab = '')
lines(1:ncol(Y), Z[2, ], col = 'black', lty = 1, lwd = 4)
abline(v = c(75, 150), h = 0, lty = 2, col = 'black', lwd = 3)




#plot(NA, xlim = c(1, ncol(Y)), ylim = c(0, 1), xlab = '', ylab = '')
plot(NA, xlim = c(1, ncol(Y)), ylim = range(Y), xlab = '', ylab = '')
col <- c('MEF' = 'blue', 'MEF_DOX' = 'blue', 'D2_EB' = 'red', 'D2_EB_DOX' = 'red')
lty <- c('MEF' = 3, 'MEF_DOX' = 1, 'D2_EB' = 3, 'D2_EB_DOX' = 1)
lapply(1:nrow(Y), function(i) lines(1:ncol(Y), Y[i, ], col = col[rownames(Y)[i]], lty = lty[rownames(Y)[i]], lwd = 4))
abline(v = c(75, 150), h = 0, lty = 2, col = 'black', lwd = 3)



#	[2017-12-01] An extended chromVAR algorithm
source('chipseq.r')




U <- matrix(sum(cvg.mnase[p]), nrow = 1, ncol = sum(n))
peaks0 <- sample.genomic.region(seqinfo(p), n = 10000, width = width(p)[1]) # sampling the background peaks
X0 <- do.call('cbind', lapply(cvg, function(x) sum(x[peaks0])))
X <- as.matrix(X %*% Diagonal(x = sum(x0) / x0))	# scale by column counts

lapply(1:nrow(d2), function(m) plot(Y[m, ], lwd = 2, type = 'l'))

cvg.mnase0 <- import(mnase.coverage.files[cell], which = reduce(peaks0), as = 'RleList')  # the depth coverage of MNase-seq


X <- (X - E) / E

X0 <- do.call('rbind', lapply(1:100, function(j) matrix(sum(cvg.mnase0[peaks0[sample(1:length(peaks), nrow(X))]]), nrow = 1, ncol = nrow(X)) %*% X))

Y0 <- do.call('rbind', lapply(1:100, function(j) matrix(sum(cvg.mnase0[peaks0[sample(1:length(peaks), nrow(X))]]), nrow = 1, ncol = nrow(X)) %*% X))



y <- U %*% X
Y0 <- do.call('rbind', lapply(1:100, function(j) matrix(sum(cvg.mnase0[peaks0[sample(1:length(peaks), nrow(X))]]), nrow = 1, ncol = nrow(X)) %*% X))



sqrt(rowSums((matrix(y, nrow = 100, ncol = length(y), byrow = TRUE) - Y0)^2))




mm <- matchMotifs(get(motif.set), peaks, genome = genome, out = 'positions')



# [2017-11-30] For each ATAC-seq file, find the motif locations
source('chipseq.r')
#for (i in 1:nrow(d)){
	cvg <- atac.coverage(d[i, 'bam.file'], genome = genome)
	peaks <- slice(cvg, lower = min.depth) 
	peaks <- as(peaks, 'GRanges')
	peaks <- reduce(peaks)	# merging the overlapping peaks
	peaks <- add.seqinfo(peaks, genome = genome)
	peaks <- trim(resize(peaks, width = atac.width, fix = 'center'))
	mm <- matchMotifs(get(motif.set), peaks, genome = genome, out = 'positions')
	mm <- add.seqinfo(mm, genome)
	cvg.mnase <- import(d[i, 'mnase.coverage.file'], which = trim(resize(peaks, width = mnase.width , fix = 'center')), as = 'RleList')  # the depth coverage of MNase-seq
	Z <- do.call('rbind', mclapply(1:length(mm), function(j){
		cat(sprintf('[%s] %d/%d processing %s\n', Sys.time(), j, length(mm), names(mm)[j]))
		if (length(mm[[j]]) > 0){
			r <- trim(resize(mm[[j]], width = mnase.width, fix = 'center'))
			MN <- as(as(cvg.mnase[r], 'RleViews'), 'matrix')
			colMeans(MN)
		}else
			rep(0, mnase.width)
	}, mc.cores = 2, mc.preschedule = FALSE))
	rownames(Z) <- names(mm)
	save(Z, file = d[i, 'nucdist.file'])
#}

#	[2017-12-01] Look the TF nucleosome distane
Z.list <- lapply(1:nrow(d), function(i){
	Z <- local({get(load(d[i, 'nucdist.file']))})
	Z <- (Z + Z[, ncol(Z):1])[, round(ncol(Z) / 2):ncol(Z)]
	Z <- t(apply(Z, 1, function(z) (z - min(z)) / (max(z) - min(z))))
	Z
})
par(mfrow = c(3, 3)); lapply(1:length(Z.list), function(i) plot(Z.list[[i]]['Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer', ]))
par(mfrow = c(3, 3)); lapply(1:length(Z.list), function(i) plot(Z.list[[i]]['ETS(ETS)/Promoter/Homer', ]))
par(mfrow = c(3, 3)); lapply(1:length(Z.list), function(i) plot(Z.list[[i]]['Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer(0.967)', ]))


#	[2017-12-01] Dynamic peaks in ATAC-seq data



mm <- matchMotifs(get(motif.set), peaks, genome = genome, out = 'positions')

X <- do.call('cbind', lapply(cvg, function(x) sum(x[peaks])))	# 



	cvg.mnase <- import(mnase.coverage.files[cell], which = trim(resize(peaks, width = mnase.width , fix = 'center')), as = 'RleList')  # the depth coverage of MNase-seq

lapply(


# [2017-11-30] First find a common peaks sets, then for each replicate, find a subset of peaks that have significant reads at these common peak sets
#cell <- 'MEF'; treatment <- 'MEF_DOX'; control <- 'MEF'
cell <- 'ESC'; treatment <- 'D2_EB_DOX'; control <- 'D2_EB'
d2 <- d[d[, 'cell'] == cell, ]

for (i in 1:length(cvg)){
	x <- sum(cvg[[i]][peaks])

	peaks2 <- reduce(as(slice(cvg[[i]], lower = min.depth), 'GRanges')





# [2017-11-30] Find the TFs that are dynamic in dynamic peaks
cell <- 'MEF'; treatment <- 'MEF_DOX'; control <- 'MEF'; peak.file <- sprintf('%s/dynamic_peaks_cell=%s_treatment=%s_control=%s.bed', project.dir, cell, treatment, control)
#cell <- 'ESC'; treatment <- 'D2_EB_DOX'; control <- 'D2_EB'; peak.file <- sprintf('%s/dynamic_peaks_cell=%s_treatment=%s_control=%s.bed', project.dir, cell, treatment, control)
d2 <- d[d[, 'cell'] == cell, ]
peak.file <- sprintf('%s/dynamic_peaks_cell=%s_treatment=%s_control=%s.RData', project.dir, cell, treatment, control)
source('chipseq.r'); peaks <- local({get(load(peak.file))})
source('chipseq.r'); cvg <- mclapply(1:nrow(d2), function(i) atac.coverage(d2[i, 'bam.file'], genome = genome), mc.cores = 4)	# coverage of ATAC-seq BAM files
dr <- !is.na(values(peaks)[['pvalue']])& values(peaks)[['pvalue']] < 0.05 & (values(peaks)[['log2FoldChange']] > log2(1.5) | values(peaks)[['log2FoldChange']] < log2(1 / 1.5))
X <- do.call('cbind', lapply(cvg, function(x) sum(x[peaks[dr]])))	# converting the Coverage object to a the read count matrix
mm <- matchMotifs(get(motif.set), peaks[dr], genome = 'mm10')	# for computeDeviations
rse <- SummarizedExperiment(assays = SimpleList(counts = X), rowRanges = peaks[dr], colData = d2[, 'group', drop = FALSE])
rse <- addGCBias(rse, genome = BSgenome.Mmusculus.UCSC.mm10)
dev <- computeDeviations(object = rse, annotations = mm)
dev.file <- sprintf('%s/computeDeviations_%s.RData', project.dir, cell)
save(dev, file = dev.file)

# [2017-11-30] For the dynamic ATAC-seq peaks, look at the MNase-seq signal surrounding each TF binding sites
#cell <- 'MEF'; treatment <- 'MEF_DOX'; control <- 'MEF'; peak.file <- sprintf('%s/dynamic_peaks_cell=%s_treatment=%s_control=%s.bed', project.dir, cell, treatment, control)
cell <- 'ESC'; treatment <- 'D2_EB_DOX'; control <- 'D2_EB'; peak.file <- sprintf('%s/dynamic_peaks_cell=%s_treatment=%s_control=%s.bed', project.dir, cell, treatment, control)
d2 <- d[d[, 'cell'] == cell, ]
peak.file <- sprintf('%s/dynamic_peaks_cell=%s_treatment=%s_control=%s.RData', project.dir, cell, treatment, control)
source('chipseq.r'); peaks <- local({get(load(peak.file))})
dr <- !is.na(values(peaks)[['pvalue']])& values(peaks)[['pvalue']] < 0.05 & values(peaks)[['log2FoldChange']] > log2(1.5)	# high in treatment
peaks2 <- peaks[dr]
mm <- matchMotifs(get(motif.set), peaks2, genome = genome, out = 'positions')
mm <- add.seqinfo(mm, genome)
cvg.mnase <- import(mnase.coverage.files[cell], which = trim(resize(peaks2, width = mnase.width , fix = 'center')), as = 'RleList')  # the depth coverage of MNase-seq
Z <- do.call('rbind', mclapply(1:length(mm), function(j){
	cat(sprintf('[%s] %d/%d processing %s\n', Sys.time(), j, length(mm), names(mm)[j]))
	r <- trim(resize(mm[[j]], width = mnase.width, fix = 'center'))
	MN <- as(as(cvg.mnase[r], 'RleViews'), 'matrix')
	colMeans(MN)
}, mc.cores = 4))
rownames(Z) <- gsub('(.+?)\\/.+', '\\1', names(mm))
Z <- (Z + Z[, ncol(Z):1])[, round(ncol(Z) / 2):ncol(Z)]
Z <- t(apply(Z, 1, function(z) (z - min(z)) / (max(z) - min(z))))
nuc.dist.file <- sprintf('%s/nuc_dist_%s.bed', project.dir, treatment)
save(Z, file = nuc.dist.file)


# []For each motif within the up-regulated ATAC-seq peaks in Dox sample, find the binding sites
cell <- 'ESC'; treatment <- 'D2_EB_DOX'; control <- 'D2_EB'; peak.file <- sprintf('%s/dynamic_peaks_cell=%s_treatment=%s_control=%s.bed', project.dir, cell, treatment, control)
d2 <- d[d[, 'cell'] == cell, ]
peak.file <- sprintf('%s/dynamic_peaks_cell=%s_treatment=%s_control=%s.RData', project.dir, cell, treatment, control)
source('chipseq.r'); peaks <- local({get(load(peak.file))})
source('chipseq.r'); cvg <- mclapply(1:nrow(d2), function(i) atac.coverage(d2[i, 'bam.file'], genome = genome), mc.cores = 4)	# coverage of ATAC-seq BAM files
up <- !is.na(values(peaks)[['pvalue']])& values(peaks)[['pvalue']] < 0.05 & values(peaks)[['log2FoldChange']] > log2(1.5)	# high in treatment




cell <- 'MEF'; treatment <- 'MEF_DOX'; control <- 'MEF'
peak.file <- sprintf('%s/dynamic_peaks_cell=%s_treatment=%s_control=%s.bed', project.dir, cell, treatment, control)
nuc.dist.file <- sprintf('%s/nuc_dist_cell=%s_treatment=%s_control=%s.bed', project.dir, cell, treatment, control)
source('chipseq.r'); Z.list <- local({get(load(nuc.dist.file))})
Z.list <- lapply(Z.list, function(Z){
	Z
})

library(pheatmap); library(gplots)
col <- colorpanel(100, low = 'blue', mid = 'black', high = 'yellow')
breaks <- seq(0, 1, length.out = 101)
pheatmap(Z, cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = FALSE, breaks = breaks, color = col, cellheight = 2, cellwidth = 2, fontsize_col = 15, fontsize_row = 9, show_colnames = FALSE, border_color = 'black', annotation_names_row = TRUE, annotation_names_col = FALSE)

col <- colorpanel(100, low = 'blue', mid = 'white', high = 'red')
dZ <- Z.list[[1]] - Z.list[[2]]
pheatmap(dZ[order(rowSums(dZ[, 1:100]), decreasing = TRUE)[1:20], ], cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = FALSE, color = col, cellheight = 9, cellwidth = 2, fontsize_col = 15, fontsize_row = 9, show_colnames = FALSE, border_color = 'black', annotation_names_row = TRUE, annotation_names_col = FALSE)
pheatmap(dZ[75:90, ], cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = FALSE, color = col, cellheight = 9, cellwidth = 2, fontsize_col = 15, fontsize_row = 9, show_colnames = FALSE, border_color = 'black', annotation_names_row = TRUE, annotation_names_col = FALSE)







# [2017-11-28] the density matrix for the dynamic NFR, using a wider range (2000 bp) for better visualization
peaks2 <- trim(resize(peaks, width = 2000, fix = 'center'))
Y <- t(sapply(Reduce('+', cvg[d[, 'group'] == 'MEF'])[peaks2[dr]], as.vector))
Y.dox <- t(sapply(Reduce('+', cvg[d[, 'group'] == 'MEF_DOX'])[peaks2[dr]], as.vector))


dhs.bed.file <- sprintf('%s/projects/sc_atac/mouse_dhs.300.bed', Sys.getenv('SHARED'))



source('chipseq.r'); cvg <- mclapply(1:nrow(d), function(i) atac.coverage(d[i, 'bam.file'], type = 'nfr', genome = 'mm10'), mc.cores = 4)	# extract the NFR reads
#source('chipseq.r'); cvg <- mclapply(1:nrow(d), function(i) atac.coverage(d[i, 'bam.file'], type = 'mononuc', genome = 'mm10'), mc.cores = 4)	# extract the NFR reads
peaks <- slice(Reduce('+', cvg), lower = 10) # regions where merged reads > 10
peaks <- reduce(as(peaks, 'GRanges'))
peaks <- add.seqinfo(peaks, genome = 'mm10')
peaks <- trim(resize(peaks, width = 500, fix = 'center'))# 500 bp regions surrounding the mononucleosomes


# --------------------------------------------------------------------------
# [2018-02-19] Look at the distribution of Etv2/CTCF binding sites over more accessible regions
# --------------------------------------------------------------------------
project <- 'Etv2_MEFs'
project.dir <- sprintf('%s/projects/%s', Sys.getenv('SHARED'), project)
dataset <- 'dataset=Naoko_version=20180220a'
source('aux.r'); se <- local({get(load(sprintf('analysis/datasets/%s.rda', dataset)))})
dds.file <- sprintf('%s/%s.dds.rda', project.dir, dataset)	# preprocessed data from DESeq2

register(MulticoreParam(4)) 
dds <- DESeqDataSet(se, design = ~ group)
dds <- DESeq(dds, parallel = TRUE)
save(dds, file = dds.file)
dds <- local({get(load(dds.file))})

library(BSgenome.Mmusculus.UCSC.mm10)
motif.set <- 'homer_pwms'
data(list = motif.set, package = 'chromVARmotifs')

dds <- addGCBias(dds, genome = BSgenome.Mmusculus.UCSC.mm10)	# add the GC content

# intervals up in the MEF OE vs MEF WT
gs.list <- list(
	c('MEF_OE', 'MEF_WT'),
	c('d2EB_OE', 'd2EB_WT'),
	c('d3EB_OE', 'd3EB_WT'),
	c('ES_OE', 'ES_WT'),
	c('Limb_OE', 'Limb_WT')
)
dev.list <- lapply(gs.list, function(gs){
	res <-	as.data.frame(results(dds, contrast = c('group', gs)))
	dr <- !is.na(res[, 'pvalue']) & res[, 'pvalue'] < 0.01 & res[, 'log2FoldChange'] > 0
	cat(sprintf('[%s] found %d differential intervals between %s and %s\n', Sys.time(), sum(dr), gs[1], gs[2]))
	m <- colData(dds)[['group']] %in% gs
	rse <- addGCBias(dds[dr, m], genome = BSgenome.Mmusculus.UCSC.mm10)
	dev <- computeDeviations(object = rse, annotations = matchMotifs(get(motif.set), rowRanges(rse), genome = 'mm10'))
	dev
})
v.list <- lapply(dev.list, function(dev) computeVariability(dev))
v2.list <- lapply(v.list, function(v) v[!grepl('ETS', rownames(v)), ])

#S <- do.call('cbind', lapply(v2.list, function(v) ranking(-log(v[, 'p_value']))))
#S <- do.call('cbind', lapply(v2.list, function(v) -log(v[, 'p_value'])))
S <- do.call('cbind', lapply(v2.list, function(v) v[, 'variability']))
rownames(S) <- v2.list[[1]][, 1]
colnames(S) <- c('MEF', 'D2EB', 'D3EB', 'ES', 'Limb')
S <- S[order(rowSums(S), decreasing = TRUE)[1:20], ]
S <- S[, c('MEF', 'ES', 'D2EB', 'D3EB', 'Limb')]

library(pheatmap); library(gplots)
col <- colorpanel(100, low = 'white', mid = 'gray', high = 'blue')
breaks <- c(-1, seq(0, 2, length.out = 99), 10)
pheatmap(S, cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = FALSE, breaks = breaks, color = col, cellheight = 15, cellwidth = 30, fontsize_col = 15, fontsize_row = 12, show_colnames = FALSE, border_color = 'black', annotation_names_row = FALSE, annotation_names_col = FALSE)

etv2 <- 'Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer(0.967)'
ctcf <- 'CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer'
ctcf2 <- 'CTCF-SatelliteElement(Zf?)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer'
gata <- 'GATA(Zf),IR4/iTreg-Gata3-ChIP-Seq(GSE20898)/Homer'
j <- 'EWS:FLI1-fusion(ETS)/SK_N_MC-EWS:FLI1-ChIP-Seq(SRA014231)/Homer'

d <- do.call('rbind', lapply(dev.list, function(dev){
	data.frame(colData(dev), t(assays(dev)[['deviations']][c(etv2, ctcf, ctcf2), ]))
}))


# --------------------------------------------------------------------------
# [2018-02-21] relationship between CTCF change on GC content
# --------------------------------------------------------------------------
project <- 'Etv2_MEFs'
project.dir <- sprintf('%s/projects/%s', Sys.getenv('SHARED'), project)
dataset <- 'dataset=Naoko_version=20180220a'
source('aux.r'); se <- local({get(load(sprintf('analysis/datasets/%s.rda', dataset)))})
dds.file <- sprintf('%s/%s.dds.rda', project.dir, dataset)	# preprocessed data from DESeq2
motif_ix.file <- sprintf('%s/%s.motif_ix.rda', project.dir, dataset)	

register(MulticoreParam(4)) 
#dds <- DESeqDataSet(se, design = ~ group)
#dds <- DESeq(dds, parallel = TRUE)
#save(dds, file = dds.file)
dds <- local({get(load(dds.file))})

library(BSgenome.Mmusculus.UCSC.mm10)
motif.set <- 'homer_pwms'
data(list = motif.set, package = 'chromVARmotifs')
dds <- addGCBias(dds, genome = BSgenome.Mmusculus.UCSC.mm10)	# add the GC content
#motif_ix <- matchMotifs(get(motif.set), rowRanges(dds), genome = 'mm10')	# add the binding information
#save(motif_ix, file = motif_ix.file)
motif_ix <- local({get(load(motif_ix.file))})


gs <-	c('d2EB_OE', 'd2EB_WT')
#gs <-	c('MEF_OE', 'MEF_WT')
res <-	as.data.frame(results(dds, contrast = c('group', gs)))
up <- !is.na(res[, 'pvalue']) & res[, 'pvalue'] < 0.01 & res[, 'log2FoldChange'] > log2(1.5)
down <- !is.na(res[, 'pvalue']) & res[, 'pvalue'] < 0.01 & res[, 'log2FoldChange'] < -log2(1.5)
ctcf <- 'CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer'
#ctcf <- 'CTCF-SatelliteElement(Zf?)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer'
etv2 <- 'Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer(0.967)'

has.ctcf <- assays(motif_ix)[['motifMatches']][, ctcf]
has.etv2 <- assays(motif_ix)[['motifMatches']][, etv2]

boxplot(list(rowData(dds[up & has.ctcf])[['bias']], rowData(dds[down & has.ctcf])[['bias']]))
boxplot(list(rowData(dds[up & has.ctcf & has.etv2])[['bias']], rowData(dds[down & has.ctcf & has.etv2])[['bias']]))

plot(rowData(dds[has.ctcf])[['bias']], res[has.ctcf, 'log2FoldChange'])


# --------------------------------------------------------------------------
# [2018-06-13] Generating bigwig files for H3K4me1, H3K27ac, H3K27me3 for ESC, MES, CP and CM
# --------------------------------------------------------------------------
project <- 'etv2_pioneer'
project.dir <- sprintf('%s/projects/%s', Sys.getenv('SHARED'), project)
dataset <- 'dataset=Wamstad_version=20180612a'
d <- read.table(sprintf('analysis/datasets/%s.tsv', dataset), sep = '\t', header = TRUE)

qvalue.cutoff <- 0.05

x <- data.frame(
	treatment = c(
		'H3K27ac_ESC', 'H3K27ac_MES', 'H3K27ac_CP', 'H3K27ac_CM',
		'H3K27me3_ESC', 'H3K27me3_MES', 'H3K27me3_CP', 'H3K27me3_CM',
		'H3K4me1_ESC', 'H3K4me1_MES', 'H3K4me1_CP', 'H3K4me1_CM',
		'H3K4me3_ESC', 'H3K4me3_MES', 'H3K4me3_CP', 'H3K4me3_CM'
	),
	control = c(
		'WCE_ESC', 'WCE_MES', 'WCE_CP', 'WCE_CM',
		'WCE_ESC', 'WCE_MES', 'WCE_CP', 'WCE_CM',
		'WCE_ESC', 'WCE_MES', 'WCE_CP', 'WCE_CM',
		'WCE_ESC', 'WCE_MES', 'WCE_CP', 'WCE_CM'
	)
)
x <- transform(x, base.name  = sprintf('%s/macs2_treatment=%s_control=%s_q=%.3e', project.dir, treatment, control, qvalue.cutoff))
library(parallel); mclapply(1:nrow(x), function(i){
	bam.files <- d[d[, 'group'] == x[i, 'treatment'], 'bam.file']
	control.files <- d[d[, 'group'] == x[i, 'control'], 'bam.file']
	source('chipseq.r'); macs2.callpeak(bam.files, base.name = x[i, 'base.name'], control.files = control.files, genome = 'mm10', broad = TRUE, nomodel = TRUE, qvalue.cutoff = qvalue.cutoff, fold.change = FALSE, update = FALSE, keep.dup = 'all', call.summits = TRUE)
}, mc.cores = 4)



# --------------------------------------------------------------------------
# [2018-06-15] Compile a list of known cardiac enhancers from two independent studeis
# (1) https://www.nature.com/articles/ng.650
# (2) dbSuper datbase
# --------------------------------------------------------------------------
project <- 'etv2_pioneer'
project.dir <- sprintf('%s/projects/%s', Sys.getenv('SHARED'), project)

source('aux.r')
library(rtracklayer)	# for liftover actions
ch <- import.chain(sprintf('%s/ucsc/mm9/mm9ToMm10.over.chain', Sys.getenv('SHARED')))	# have to be uncompressed file

# https://www.ncbi.nlm.nih.gov/pubmed/20729851
x <- read.table(sprintf('%s/heart_enhancers.tsv', project.dir), header = TRUE, sep = '\t')	# the coordinates are in mm9
chrom <- gsub('^(.+) (.+) (.+)$', '\\1', x[, 1])
from <- as.numeric(gsub('^(.+) (.+) (.+)$', '\\2', x[, 1]))
to <- as.numeric(gsub('^(.+) (.+) (.+)$', '\\3', x[, 1]))
x <- GRanges(seqnames = chrom, range = IRanges(from, to))
x <- add.seqinfo(x, genome = 'mm9')
seqlevelsStyle(x) <- 'UCSC'  # necessary
x <- unlist(liftOver(x, ch))
x <- add.seqinfo(x, genome = 'mm10')

x2 <- read.table(sprintf('%s/20180618_Heart.bed', project.dir), header = FALSE, sep = '\t')
x2 <- GRanges(seqnames = x2[, 1], range = IRanges(x2[, 2], x2[, 3]))
x2 <- add.seqinfo(x2, genome = 'mm9')
seqlevelsStyle(x2) <- 'UCSC'  # necessary
x2 <- unlist(liftOver(x2, ch))
x2 <- add.seqinfo(x2, genome = 'mm10')

x3 <- read.table(sprintf('%s/20180618_E14.5_Heart.bed', project.dir), header = FALSE, sep = '\t')
x3 <- GRanges(seqnames = x3[, 1], range = IRanges(x3[, 2], x3[, 3]))
x3 <- add.seqinfo(x3, genome = 'mm9')
seqlevelsStyle(x3) <- 'UCSC'  # necessary
x3 <- unlist(liftOver(x3, ch))
x3 <- add.seqinfo(x3, genome = 'mm10')

enhancer <- reduce(c(x, x2, x3))	# merged set of mouse heart enhancers from two studies
values(enhancer)[['source']] <- cbind(
	blow = 1:length(enhancer) %in% as.matrix(findOverlaps(enhancer, x))[, 1], 
	dbsuper_adult_heart = 1:length(enhancer) %in% as.matrix(findOverlaps(enhancer, x2))[, 1],
	dbsuper_E14_heart = 1:length(enhancer) %in% as.matrix(findOverlaps(enhancer, x3))[, 1]
)
# adding the conservation scores
phastcons.file <- sprintf('%s/igenome/Mus_musculus/mm10.60way.phastCons.bw', Sys.getenv('SHARED'))
cvg <- import(phastcons.file, which = enhancer, as = 'RleList')	# the conservation coverage
phastCons <- mean(cvg[enhancer])
values(enhancer)[['phastCons']] <- phastCons

dataset <- 'dataset=cardiac_enhancers_version=20180619a'
save(enhancer, file = sprintf('analysis/datasets/%s.rda', dataset))


# --------------------------------------------------------------------------
# [2018-06-19] Generate object from ChIP-Atlas database (mouse)
# The Gata2/Gata4 ChIP-seq data were downloaded from http://chip-atlas.org/peak_browser
# copying BED files downloaded from ChIP-Atlas to MSI (e.g.):
# scp ~/Desktop/Oth.ALL.05.Tal1.AllCell.bed gongx030@login.msi.umn.edu:/home/garrydj/gongx030/rlib/analysis/datasets/chip_atlas/mm9
# --------------------------------------------------------------------------
source('aux.r')
library(rtracklayer)	# for liftover actions
# https://master.bioconductor.org/packages/release/workflows/vignettes/liftOver/inst/doc/liftov.html
ch <- import.chain(sprintf('%s/ucsc/mm9/mm9ToMm10.over.chain', Sys.getenv('SHARED')))	# have to be uncompressed file
chip.atlas.dir <- sprintf('analysis/datasets/chip_atlas/mm9', Sys.getenv('SHARED'))

#bed.file <- sprintf('%s/Oth.ALL.05.Gata2.AllCell.bed', chip.atlas.dir); antibody <- 'Gata2'
#bed.file <- sprintf('%s/Oth.ALL.05.Gata4.AllCell.bed', chip.atlas.dir); antibody <- 'Gata4'
#bed.file <- sprintf('%s/Oth.ALL.05.Tal1.AllCell.bed', chip.atlas.dir); antibody <- 'Tal1'

x <- read.table(bed.file, sep = '\t', header = FALSE, quote = '', comment.char = '', skip = 1)
cell_group <- gsub('.+Cell%20group=(.+?);.+', '\\1', x[, 4]); cell_group <- gsub('%20', '_', cell_group)
id <- gsub('ID=(.+?);.+', '\\1', x[, 4])
name <- gsub('.+;Name=(.+?);.+', '\\1', x[, 4]); name <- gsub('%20', '_', name)
title <- gsub('.+;Title=(.+?);.+', '\\1', x[, 4]); title <- gsub('%20', '_', title); title <- gsub('%3B', '', title)
source_name <- gsub('.+source_name=(.+?);.+', '\\1', x[, 4]); source_name <- gsub('%20', '_', source_name)
gr <- GRanges(seqnames = x[, 1], range = IRanges(x[, 2], x[, 3]), antibody = antibody, cell_group = cell_group, id = id, name = name, title = title, source_name = source_name)
gr <- add.seqinfo(gr, genome = 'mm9')
seqlevelsStyle(gr) <- 'UCSC'  # necessary
gr <- unlist(liftOver(gr, ch))
gr <- add.seqinfo(gr, genome = 'mm10')
save(gr, file = sprintf('analysis/datasets/chip_atlas/mm10/%s.rda', antibody))



# --------------------------------------------------------------------------
# [2018-06-15] Compile a list of known cardiac enhancers from two independent studeis
# --------------------------------------------------------------------------
values(enhancer)$has_gata2 = 	1:length(enhancer) %in% as.matrix(findOverlaps(enhancer, gata2))[, 1]
values(enhancer)$has_gata4 = 	1:length(enhancer) %in% as.matrix(findOverlaps(enhancer, gata4))[, 1]
values(enhancer)$has_etv2 = 	1:length(enhancer) %in% as.matrix(findOverlaps(enhancer, etv2))[, 1]

# find known cardiac genes and associated TSS
library(org.Mm.eg.db)
sym <- unique(unlist(as.list(org.Mm.egSYMBOL)))
source('go.r'); G <- go.matrix(sym, filter = 'external_gene_name', domain = 'BP', method = 'ncbi')
heart.genes <- names(which(G[, 'GO:0007507']))
library(biomaRt)
mart <- useMart(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl', host = 'www.ensembl.org')
bm <- getBM(attributes = c('external_gene_name', 'ensembl_transcript_id', 'chromosome_name', 'transcription_start_site', 'strand'), filters = 'external_gene_name', values = heart.genes, mart = mart)

# add the expression change on Etv2 or Gata2 induction to these heart genes
z <- read.table(sprintf('%s/Etv2_Gata2_RNAseq.tsv', project.dir), header = TRUE, sep = '\t')
z <- data.frame(symbol = sprintf('%s%s', substr(z[, 1], 1, 1), tolower(substr(z[, 1], 2, nchar(z[, 1])))), WT = z[, 2], Etv2_OE = z[, 4], Gata2_OE = z[, 5], FC_Etv2_OE = log2((z[, 4] + 1) / (z[, 2] + 1)), FC_Gata2_OE = log2((z[, 5] + 1) / (z[, 2] + 1)))
z <- z[z[, 'Etv2_OE'] > 5 | z[, 'WT'] > 5 | z[, 'Gata2_OE'] > 5, ]
z <- merge(bm, z, by.x = 'external_gene_name', by.y = 'symbol')

heart.tx <- GRanges(
	seqnames = sprintf('chr%s', z[, 'chromosome_name']), range = IRanges(z[, 'transcription_start_site'], z[, 'transcription_start_site']), strand = z[, 'strand'], 
	ensembl_transcript_id = z[, 'ensembl_transcript_id'], 
	external_gene_name = z[, 'external_gene_name'],
	WT = z[, 'WT'],
	Etv2_OE = z[, 'Etv2_OE'],
	Gata2_OE = z[, 'Gata2_OE'],
	FC_Etv2_OE = z[, 'FC_Etv2_OE'],
	FC_Gata2_OE = z[, 'FC_Gata2_OE']
)
heart.tx <- add.seqinfo(heart.tx, genome = 'mm10')
mm <- as.matrix(findOverlaps(trim(resize(heart.tx, width = 1000000, fix = 'center')), enhancer))

xx <- cbind(as.data.frame(enhancer[mm[, 2]]), as.data.frame(heart.tx[mm[, 1]]), distance = distance(enhancer[mm[, 2]], heart.tx[mm[, 1]]))
etv2 <- read.table(sprintf('%s/Etv2_chipseq.bed', project.dir), sep = '\t', header = FALSE)
etv2 <- GRanges(seqnames = etv2[, 1], range = IRanges(etv2[, 2], etv2[, 3]))
mm <- as.matrix(findOverlaps(etv2, GRanges(seqnames = xx[, 1], range = IRanges(xx[, 2], xx[, 3]))))
xx <- cbind(xx[mm[, 2], ], as.data.frame(etv2[mm[, 1]]))
xx <- xx[xx[, 'FC_Etv2_OE'] < -0.25, ]

write.table(xx, sprintf('%s/heart_enancers_Etv2_Gata2.tsv', project.dir), sep = '\t', quote = FALSE, row.names = FALSE)

gplots::venn(list(Gata2 = which(values(enhancer)$has_gata2), Gata4 = which(values(enhancer)$has_gata4), Etv2 = which(values(enhancer)$has_etv2)))


# --------------------------------------------------------------------------
# [2018-06-17] Competition of Gata2 and Gata4 binding sites
# The Gata2/Gata4 ChIP-seq data were downloaded from http://chip-atlas.org/peak_browser
# --------------------------------------------------------------------------
gata2 <- read.table(sprintf('%s/Oth.ALL.05.Gata2.AllCell.bed', project.dir), sep = '\t', header = FALSE, quote = '', comment.char = '')
gata2 <- gata2[gsub('.+Cell%20group=(.+?);.+', '\\1', gata2[, 4]) %in% c('Blood', 'Pluripotent%20stem%20cell'), ]	# mesoderm lineage
gata2 <- GRanges(seqnames = gata2[, 1], range = IRanges(gata2[, 2], gata2[, 3]))
gatq2 <- add.seqinfo(gata2, genome = 'mm9')
seqlevelsStyle(gata2) <- 'UCSC'  # necessary
gata2 <- unlist(liftOver(gata2, ch))
gata2 <- add.seqinfo(gata2, genome = 'mm10')
gata2 <- reduce(gata2)

gata4 <- read.table(sprintf('%s/Oth.ALL.05.Gata4.AllCell.bed', project.dir), sep = '\t', header = FALSE, quote = '', comment.char = '')
gata4 <- gata4[gsub('.+Cell%20group=(.+?);.+', '\\1', gata4[, 4]) %in% c('Cardiovascular', 'Pluripotent%20stem%20cell'), ]	# mesoderm lineage
gata4 <- GRanges(seqnames = gata4[, 1], range = IRanges(gata4[, 2], gata4[, 3]))
gata4 <- add.seqinfo(gata4, genome = 'mm9')
seqlevelsStyle(gata4) <- 'UCSC'  # necessary
gata4 <- unlist(liftOver(gata4, ch))
gata4 <- add.seqinfo(gata4, genome = 'mm10')
gata4 <- reduce(gata4)

gata <- reduce(c(gata2, gata4))

xx <- cbind(
	Gata2 = 1:length(gata) %in% as.matrix(findOverlaps(gata, gata2))[, 1],
	Gata4 = 1:length(gata) %in% as.matrix(findOverlaps(gata, gata4))[, 1],
	Gata2_Gata4 = 1:length(gata) %in% as.matrix(findOverlaps(gata, gata4))[, 1] & 1:length(gata) %in% as.matrix(findOverlaps(gata, gata2))[, 1],
	Gata2_enhancer = 1:length(gata) %in% as.matrix(findOverlaps(gata, gata2))[, 1] & 1:length(gata) %in% as.matrix(findOverlaps(gata, enhancer))[, 1],
	Gata4_enhancer = 1:length(gata) %in% as.matrix(findOverlaps(gata, gata4))[, 1] & 1:length(gata) %in% as.matrix(findOverlaps(gata, enhancer))[, 1],
	Gata2_Gata4_enhancer = 1:length(gata) %in% as.matrix(findOverlaps(gata, gata4))[, 1] & 1:length(gata) %in% as.matrix(findOverlaps(gata, gata2))[, 1] & 1:length(gata) %in% as.matrix(findOverlaps(gata, enhancer))[, 1]
)
y <- colSums(xx)



# --------------------------------------------------------------------------
# [2019-03-01] Read ATAC-seq BAM files and convert into intervals of reads
# --------------------------------------------------------------------------
dataset <- 'dataset=Etv2ATAC_version=20190228a'
genome <- 'mm10'; atac.width <- 500; min.depth <- 20
source('sra.r'); d <- read.dataset(dataset)
d <- transform(d, sra.run.dir = sra.run.dir(run), bam.file = sprintf('%s/%s.dedup.bam', sra.run.result.dir(run), run))
source('chipseq.r'); se <- atac.read.peaks(d[, 'bam.file'], genome = genome, min.depth = min.depth, atac.width = atac.width, mc.cores = 2)
colData(se) <- DataFrame(d)
peak_file <- sprintf('analysis/etv2_pioneer/data/etv2_atacseq_peaks.rds')
saveRDS(se, file = peak_file)


# --------------------------------------------------------------------------
# [2019-05-06] Read the sample level ATAC-seq data into a count matrix
# --------------------------------------------------------------------------
library(chromVAR)
library(BiocParallel)
register(MulticoreParam(2)) # Use 8 cores
peak_file <- sprintf('analysis/etv2_pioneer/data/etv2_atacseq_peaks.rds')
peaks <- readRDS(peak_file)
source('analysis/etv2_pioneer/helper.r'); d <- read_ATAC_dataset(sample_level = FALSE)
source('chipseq.r'); se <- get_counts(d[, 'bam.file'], peaks)
colData(se)$group <- d$group

library(SummarizedExperiment)
library(BSgenome.Mmusculus.UCSC.mm10)
library(chromVARmotifs) # https://github.com/GreenleafLab/chromVARmotifs
library(motifmatchr)

motif.set <- 'homer_pwms'
data(list = motif.set, package = 'chromVARmotifs')
se <- addGCBias(se, genome = BSgenome.Mmusculus.UCSC.mm10)
motif_ix <- matchMotifs(get(motif.set), se, genome = 'mm10')
dev <- computeDeviations(object = se, annotations = motif_ix)
v <- computeVariability(dev)
												                                  
h <- v[, 'p_value_adj'] < 0.05
Z <- assays(dev)$z[h, ]
pca <- prcomp(Z)
source('analysis/etv2_pioneer/helper.r'); bg <- atacseq_group2bg(colData(se)$group)
source('analysis/etv2_pioneer/helper.r'); pch <- atacseq_group2pch(colData(se)$group)
source('analysis/etv2_pioneer/helper.r'); col <- atacseq_group2col(colData(se)$group)
plot(pca$rotation[, 1], pca$rotation[, 2], pch = pch, bg = bg, col = col, xaxt = 'n', yaxt = 'n', xlab = sprintf('PC1(%.1f%%)',  100 * summary(pca)$importance[2, 1]), ylab = sprintf('PC2(%.1f%%)',  100 * summary(pca)$importance[2, 2]), cex = 2, lwd = 2)
plot(pca$rotation[, 3], pca$rotation[, 2], pch = pch, bg = bg, col = col, xaxt = 'n', yaxt = 'n', xlab = sprintf('PC3(%.1f%%)',  100 * summary(pca)$importance[2, 3]), ylab = sprintf('PC2(%.1f%%)',  100 * summary(pca)$importance[2, 2]), cex = 2, lwd = 2)


# --------------------------------------------------------------------------
# [2019-05-06] Read the condition level ATAC-seq data into a count matrix
# And do the chromVAR analysis
# --------------------------------------------------------------------------
library(chromVAR)
library(BiocParallel)
register(MulticoreParam(4)) # Use 8 cores
peak_file <- sprintf('analysis/etv2_pioneer/data/etv2_atacseq_peaks.rds')
peaks <- readRDS(peak_file)
source('chipseq.r'); se <- get_counts(d2[, 'bam_file'], peaks)


# --------------------------------------------------------------------------
# [2019-03-01] chromVAR analysis of Etv2 ATAC-seq in MEFs and ES/EB
# --------------------------------------------------------------------------
library(pheatmap)
h <- unique(unlist(lapply(1:ncol(Z), function(m) order(Z[, m], decreasing = TRUE)[1:30])))
Z2 <- Z[h, ]
rownames(Z2) <- gsub('(.+?)/.+', '\\1', rownames(Z2)); colnames(Z2) <- sprintf('S:%d', 1:ncol(Z2))
df <- data.frame(group = colData(se)$group)
rownames(df) <- colnames(Z2)
collist <- list(group = atacseq_group2bg(colData(se)$group))
col <- colorpanel(100, low = 'blue', mid = 'white', high = 'red')
breaks <- c(-100, seq(-20, 20, length.out = 99), 100)
pheatmap(Z2, cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = TRUE, breaks = breaks, fontsize_row = 9, annotation_col = df, color = col, cellheight = 9, cellwidth = 15, show_colnames = FALSE, border_color = 'black', annotation_colors = collist, annotation_names_row = FALSE, annotation_names_col = FALSE)



# --------------------------------------------------------------------------
# [2019-03-08] calling Etv2 peaks in MEF reprogramming
# --------------------------------------------------------------------------
dataset <- 'dataset=Etv2ChIPseq_version=20190307a'
source('sra.r'); d <- read.dataset(dataset, touch = FALSE)
treatment <- rbind(
	d[, 'group'] %in% c('MEF_Dox_D1_Etv2_rep1', 'MEF_Dox_D1_Etv2_rep2'),
	d[, 'group'] %in% c('MEF_Dox_D2_Etv2_rep1', 'MEF_Dox_D2_Etv2_rep2'),
	d[, 'group'] %in% c('MEF_Dox_D7_Etv2_rep1', 'MEF_Dox_D7_Etv2_rep2'),
	d[, 'group'] %in% c('Limb_E10.5_GFPpos_Etv2'),
	d[, 'group'] %in% c('Limb_E10.5_GFPneg_Etv2'),

	d[, 'group'] %in% c('MEF_NoDox_Brg1_rep1', 'MEF_NoDox_Brg1_rep2'),
	d[, 'group'] %in% c('MEF_Dox_D1_Brg1_rep1', 'MEF_Dox_D1_Brg1_rep2'),
	d[, 'group'] %in% c('MEF_Dox_D2_Brg1_rep1', 'MEF_Dox_D2_Brg1_rep2'),
	d[, 'group'] %in% c('MEF_Dox_D7_Brg1_rep1', 'MEF_Dox_D7_Brg1_rep2'),
	d[, 'group'] %in% c('Limb_E10.5_GFPpos_Brg1'),
	d[, 'group'] %in% c('Limb_E10.5_GFPneg_Brg1'),

	d[, 'group'] %in% c('MEF_NoDox_H3K27ac_rep1', 'MEF_NoDox_H3K27ac_rep2'),
	d[, 'group'] %in% c('MEF_Dox_D1_H3K27ac_rep1', 'MEF_Dox_D1_H3K27ac_rep2'),
	d[, 'group'] %in% c('MEF_Dox_D2_H3K27ac_rep1', 'MEF_Dox_D2_H3K27ac_rep2'),
	d[, 'group'] %in% c('MEF_Dox_D7_H3K27ac_rep1', 'MEF_Dox_D7_H3K27ac_rep2')
)

control <- rbind(
	d[, 'group'] %in% c('MEF_Dox_D1_input_rep1'),
	d[, 'group'] %in% c('MEF_Dox_D2_input_rep1'),
	d[, 'group'] %in% c('MEF_Dox_D7_input_rep1'),
	d[, 'group'] %in% c('Limb_E10.5_GFPpos_input'),
	d[, 'group'] %in% c('Limb_E10.5_GFPneg_input'),

	d[, 'group'] %in% c('MEF_NoDox_input_rep1'),
	d[, 'group'] %in% c('MEF_Dox_D1_input_rep1'),
	d[, 'group'] %in% c('MEF_Dox_D2_input_rep1'),
	d[, 'group'] %in% c('MEF_Dox_D7_input_rep1'),
	d[, 'group'] %in% c('Limb_E10.5_GFPpos_input'),
	d[, 'group'] %in% c('Limb_E10.5_GFPneg_input'),

	d[, 'group'] %in% c('MEF_NoDox_input_rep1'),
	d[, 'group'] %in% c('MEF_Dox_D1_input_rep1'),
	d[, 'group'] %in% c('MEF_Dox_D2_input_rep1'),
	d[, 'group'] %in% c('MEF_Dox_D7_input_rep1')
)

base.name <- sprintf('analysis/etv2_pioneer/results/%s', c(
		'MEF_Dox_D1_Etv2', 
		'MEF_Dox_D2_Etv2',
		'MEF_Dox_D7_Etv2',
		'Limb_E10.5_GFPpos_Etv2',
		'Limb_E10.5_GFPneg_Etv2',
		'MEF_NoDox_Brg1',
		'MEF_Dox_D1_Brg1',
		'MEF_Dox_D2_Brg1',
		'MEF_Dox_D7_Brg1',
		'Limb_E10.5_GFPpos_Brg1',
		'Limb_E10.5_GFPneg_Brg1',
		'MEF_NoDox_H3K27ac',
		'MEF_Dox_D1_H3K27ac',
		'MEF_Dox_D2_H3K27ac',
		'MEF_Dox_D7_H3K27ac'
))
broad <- rep(FALSE, length(base.name))
broad[grep('H3K27ac', base.name)] <- TRUE
call_submit <- rep(TRUE, length(base.name))
call_submit[grep('H3K27ac', base.name)] <- FALSE 
pileup_file <- sprintf('%s_treat_pileup.bw', base.name)
s3_pileup_file <- gsub('.+/(.+)', 's3://etv2_pioneer/\\1', pileup_file)
s3_pileup_public_file <- gsub('.+/(.+)', 'https://s3.msi.umn.edu/etv2_pioneer/\\1', pileup_file)
y <- DataFrame(treatment = I(treatment), control = I(control), base.name = base.name, pileup_file = pileup_file, s3_pileup_file = s3_pileup_file, s3_pileup_public_file = s3_pileup_public_file)

treatment_files <- lapply(1:nrow(y), function(i) d[y[i, 'treatment'], 'bam.file'])
# all(unlist(lapply(treatment_files, file.exists)))
control_files <- lapply(1:nrow(y), function(i) d[y[i, 'control'], 'bam.file'])
# all(unlist(lapply(control_files, file.exists)))


source('chipseq.r'); 
lapply(which(!call_submit), function(i){
	macs2.callpeak(treatment_files[[i]], y[i, 'base.name'], control_files[[i]], genome = 'mm10', broad = broad[i], qvalue.cutoff = 0.05, fold.change = FALSE, update = TRUE, call.summits = call_submit[i])
})

source('s3.r'); 
res <- lapply(1:nrow(y), function(i){
	s3.backup(y[i, 'pileup_file'], y[i, 's3_pileup_file'], make_public = TRUE)
})


# --------------------------------------------------------------------------
# [2019-03-08] Find the nucleosome positions of Etv2 ChIP-seq peaks
# on MEF Etv2 Dox samples
# [2019-04-18] Line plot surrounding the Etv2 ChIP-seq peaks on D1
# --------------------------------------------------------------------------
source('chipseq.r')
peaks_file <- sprintf('analysis/etv2_pioneer/results/MEF_Dox_D1_Etv2_summits.bed')
gr <- read.table(peaks_file, header = FALSE, sep = '\t')
gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]), peak_id = gr[, 4], log10pvalue = -gr[, 5])
gr <- add.seqinfo(gr, genome = 'mm10')
gr <- gr[values(gr)[['log10pvalue']] < -10]
gr <- resize(gr, width = 2000, fix = 'center')	# extending centered at the summit

# download the MEF nucleosome pileup bw
bw_files <- c(
#	D1_Brg1 = 'analysis/etv2_pioneer/results/MEF_Dox_D1_Brg1_treat_pileup.bw',
#	NoDox_Brg1 = 'analysis/etv2_pioneer/results/MEF_NoDox_Brg1_treat_pileup.bw',
	NoDox_Brg1 = 'analysis/etv2_pioneer/results/MEF_Brg1_GSE90893_treat_pileup.bw',
	D1_H3K27ac = 'analysis/etv2_pioneer/results/MEF_Dox_D1_H3K27ac_treat_pileup.bw',
	NoDox_H3K27ac = 'analysis/etv2_pioneer/results/MEF_NoDox_H3K27ac_treat_pileup.bw',
	D1_Etv2 = 'analysis/etv2_pioneer/results/MEF_Dox_D1_Etv2_treat_pileup.bw'
)

Z <- do.call('rbind', lapply(bw_files, function(bw_file){
	cat(sprintf('reading %s\n', bw_file))
	cvg <- import(bw_file , which = trim(reduce(gr)), as = 'RleList')  
	Z <- as(as(cvg[gr], 'RleViews'), 'matrix')
	Z <- t(apply(Z, 1, function(x) (x - min(x)) / (max(x) - min(x))))
	Z[is.na(Z)] <- 0
	colMeans(Z)
}))
Z2 <- t(apply(Z, 1, function(x) (x - min(x))))

plot(Z2['D1_Etv2', ], type = 'l', lwd = 2, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
#lines(Z2['D1_Brg1', ], lwd = 2, col = 'red')
lines(Z2['NoDox_Brg1', ], lwd = 2, lty = 2, col = 'red')
lines(Z2['D1_H3K27ac', ], lwd = 2, col = 'blue')
lines(Z2['NoDox_H3K27ac', ], lwd = 2, lty = 2, col = 'blue')



library(TxDb.Mmusculus.UCSC.mm10.knownGene)
pmt <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene)
is_promoter <- 1:length(gr) %in% as.matrix(findOverlaps(gr, pmt))[, 1]

H3K27ac_peak_file <- 'analysis/etv2_pioneer/results/MEF_NoDox_H3K27ac_peaks.broadPeak'
gr_H3K27ac <- read.table(H3K27ac_peak_file, header = FALSE, sep = '\t')
gr_H3K27ac <- GRanges(seqnames = gr_H3K27ac[, 1], range = IRanges(gr_H3K27ac[, 2], gr_H3K27ac[, 3]), peak_id = gr_H3K27ac[, 4], log10pvalue = -gr_H3K27ac[, 5])
gr_H3K27ac <- add.seqinfo(gr_H3K27ac, genome = 'mm10')
gr_H3K27ac <- resize(gr_H3K27ac, width = 5000, fix = 'center')	# extending centered at the summit
is_enhancer <- 1:length(gr) %in% as.matrix(findOverlaps(gr, gr_H3K27ac))[, 1]

i <- !is_promoter & !is_enhancer	# non promoter regions and not marked by H3K27ac at MEF


# --------------------------------------------------------------------------
# [2019-03-08] Annotate ChIP-seq peaks
# --------------------------------------------------------------------------
peaks_files <- sprintf('analysis/etv2_pioneer/results/%s_peaks.narrowPeak', c(
		'MEF_Dox_D1_Etv2', 
		'MEF_Dox_D2_Etv2',
		'MEF_Dox_D7_Etv2',
		'Limb_E10.5_GFPpos_Etv2',
		'Limb_E10.5_GFPneg_Etv2',
		'MEF_NoDox_Brg1',
		'MEF_Dox_D1_Brg1',
		'MEF_Dox_D2_Brg1',
		'MEF_Dox_D7_Brg1',
		'Limb_E10.5_GFPpos_Brg1',
		'Limb_E10.5_GFPneg_Brg1',
		'MEF_NoDox_H3K27ac',
		'MEF_Dox_D1_H3K27ac',
		'MEF_Dox_D2_H3K27ac',
		'MEF_Dox_D7_H3K27ac'
))
source('chipseq.r'); homer(peaks_files, genome = 'mm10', mc.cores = 8)
output_dirs <- sprintf('%s_homer', peaks_files)
source('s3.r'); s3.backup(output_dirs, 's3://etv2_pioneer/homer/', make_public = TRUE)


# --------------------------------------------------------------------------
# [2019-03-12] Look at the % of overlap of Etv2 ChIP-seq peaks between
# D1, D3 and D7
# [2019-03-14] A unified Etv2 summits across D1, D2 and D7
# --------------------------------------------------------------------------
library(SummarizedExperiment)
source('chipseq.r')
peaks_files <- sprintf('analysis/etv2_pioneer/results/%s_summits.bed', c(
		'MEF_Dox_D1_Etv2', 
		'MEF_Dox_D2_Etv2',
		'MEF_Dox_D7_Etv2'
))
grl <- lapply(peaks_files, function(file){
	gr <- read.table(file, header = FALSE, sep = '\t')
	gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]), peak_id = gr[, 4], log10pvalue = -gr[, 5])
	gr <- resize(gr, width = 250, fix = 'center')
	gr <- gr[values(gr)[['log10pvalue']] < -10]
	gr <- add.seqinfo(gr, genome = 'mm10')
	gr
})
gr0 <- Reduce('union', grl)
for (j in 1:length(grl)){
	mm <- as.matrix(findOverlaps(gr0, grl[[j]]))
	values(grl[[j]])$peak_group[mm[, 2]] <- mm[, 1] 
}
gr <- Reduce('c', grl)
sp <- split(1:length(gr), list(factor(values(gr)$peak_group)))
log10pvalue <- values(gr)$log10pvalue
keep <- sapply(sp, function(i) i[which.min(log10pvalue[i])])
gr <- gr[keep]
Y <- do.call('cbind', lapply(grl, function(x) 1:length(gr) %in% as.matrix(findOverlaps(gr, x))[, 1]))
values(gr)$present <- Y
file <- sprintf('analysis/etv2_pioneer/results/MEF_Etv2_peaks.rds')
saveRDS(gr, file)


# --------------------------------------------------------------------------
# [2019-03-12] Look at the % of overlap of H3K27ac ChIP-seq peaks between
# D1, D3 and D7
# --------------------------------------------------------------------------
source('chipseq.r')
peaks_files <- sprintf('analysis/etv2_pioneer/results/%s_peaks.broadPeak', c(
	'MEF_NoDox_H3K27ac',
	'MEF_Dox_D1_H3K27ac',
	'MEF_Dox_D2_H3K27ac',
	'MEF_Dox_D7_H3K27ac'
))
grl <- lapply(peaks_files, function(file){
	gr <- read.table(file, header = FALSE, sep = '\t')
	gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]), peak_id = gr[, 4], log10pvalue = -gr[, 5])
	gr <- resize(gr, width = 500, fix = 'center')
	gr <- add.seqinfo(gr, genome = 'mm10')
	gr
})
gr0 <- Reduce('union', grl)
Y <- do.call('cbind', lapply(1:length(grl), function(j){
	 1:length(gr0) %in% as.matrix(findOverlaps(gr0, grl[[j]]))
}))
class(Y) <- 'numeric'

library(pheatmap); library(gplots)
code <- sprintf('%d%d%d%d', Y[, 1], Y[, 2], Y[, 3], Y[, 4])
table(code)


# --------------------------------------------------------------------------
# [2019-03-12] Calling peaks of ATAC-seq data
# [2019-05-16] These are the sample level data; 
# !!! there might be some wrong with these codes, since the ATAC-seq in MEF appears not changing !!!
# --------------------------------------------------------------------------
dataset <- 'dataset=Etv2ATAC_version=20190228a'
source('sra.r'); d <- read.dataset(dataset, touch = FALSE)
treatment <- rbind(
	d[, 'group'] %in% c('MEF_NoDox_rep1', 'MEF_NoDox_rep2'),
	d[, 'group'] %in% c('MEF_Dox_D1_rep1', 'MEF_Dox_D1_rep2'),
	d[, 'group'] %in% c('MEF_Dox_D2_rep1', 'MEF_Dox_D2_rep2'),
	d[, 'group'] %in% c('MEF_Dox_D7_rep1', 'MEF_Dox_D7_rep2'),
	d[, 'group'] %in% c('MEF_Dox_D7_Flk1pos_rep1', 'MEF_Dox_D7_Flk1pos_rep2'),
	d[, 'group'] %in% c('EB_Dox_D25_rep1', 'EB_Dox_D25_rep2'),
	d[, 'group'] %in% c('EB_Dox_D25_Flk1pos_rep1', 'MEF_Dox_D7_Flk1pos_rep2'),
	d[, 'group'] %in% c('EB_NoDox_D25_rep1', 'EB_NoDox_D25_rep2')
)
base.name <- sprintf('analysis/etv2_pioneer/results/ATAC_%s', c(
		'MEF_NoDox', 
		'MEF_Dox_D1', 
		'MEF_Dox_D2', 
		'MEF_Dox_D7', 
		'MEF_Dox_D7_Flk1pos', 
		'EB_Dox_D25', 
		'EB_Dox_D25_Flk1pos', 
		'EB_NoDox_D25'
))
pileup_file <- sprintf('%s_treat_pileup.bw', base.name)
s3_pileup_file <- gsub('.+/(.+)', 's3://etv2_pioneer/\\1', pileup_file)
s3_pileup_public_file <- gsub('.+/(.+)', 'https://s3.msi.umn.edu/etv2_pioneer/\\1', pileup_file)
y <- DataFrame(treatment = I(treatment), base.name = base.name, pileup_file = pileup_file, s3_pileup_file = s3_pileup_file, s3_pileup_public_file = s3_pileup_public_file)

treatment_files <- lapply(1:nrow(y), function(i) d[y[i, 'treatment'], 'bam.file'])
all(unlist(lapply(treatment_files, file.exists)))

source('chipseq.r'); 
lapply(1:nrow(y), function(i){
	macs2.callpeak(treatment_files[[i]], y[i, 'base.name'], format = 'BAMPE', genome = 'mm10', broad = FALSE, qvalue.cutoff = 0.05, fold.change = FALSE, update = TRUE, call.summits = TRUE, shift = -100, extsize = 200)
})

source('s3.r'); 
res <- lapply(1:nrow(y), function(i){
	s3.backup(y[i, 'pileup_file'], y[i, 's3_pileup_file'], make_public = TRUE)
})


# --------------------------------------------------------------------------
# [2019-03-13] Look at the Etv2 motifs of each Etv2 binding site cluster
# --------------------------------------------------------------------------
library(SummarizedExperiment)
source('chipseq.r')
file <- sprintf('analysis/etv2_pioneer/results/Etv2_peaks_cluster.rds')
gr <- readRDS(file)
x <- as.data.frame(gr[values(gr)$cluster == '100'])
peak_file <- sprintf('analysis/etv2_pioneer/results/Etv2_peaks_cluster=100.bed')
write.table(x, peak_file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
source('chipseq.r'); homer(peak_file, genome = 'mm10', mc.cores = 8)
output_dir <- sprintf('%s_homer', peak_file)
source('s3.r'); s3.backup(output_dir, 's3://etv2_pioneer/results/', make_public = TRUE)



# --------------------------------------------------------------------------
# [2019-03-13] Find the nuleosome free and mononucleosome reads from ATAC-seq data
# at NoDox and D1 MEF, centered at the Etv2 binding sites at D1
# [2019-03-14] NucleoATAC identified too few NFR regions
# --------------------------------------------------------------------------
source('chipseq.r')
peak_file <- sprintf('analysis/etv2_pioneer/results/%s_summits.bed', 'MEF_Dox_D1_Etv2')
gr <- read.table(peak_file, header = FALSE, sep = '\t')
gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]), peak_id = gr[, 4], log10pvalue = -gr[, 5])
gr <- gr[values(gr)[['log10pvalue']] < -10]
gr <- resize(gr, width = 2000, fix = 'center')	# extending centered at the summit
gr <- add.seqinfo(gr, genome = 'mm10')
gr <- reduce(trim(gr))
peak_file2 <- sprintf('analysis/etv2_pioneer/results/%s_summits_w2000.bed', 'MEF_Dox_D1_Etv2')
write.table(as.data.frame(gr)[, 1:3], peak_file2, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

source('chipseq.r')
peak_file2 <- sprintf('analysis/etv2_pioneer/results/%s_summits_w2000.bed', 'MEF_Dox_D1_Etv2')
dataset <- 'dataset=Etv2ATAC_version=20190228a'
source('sra.r'); d <- read.dataset(dataset, touch = FALSE)
#groups <- c('MEF_NoDox_rep1', 'MEF_NoDox_rep2'); output_prefix <- 'ATAC_MEF_NoDox'
#groups <- c('MEF_Dox_D1_rep1', 'MEF_Dox_D1_rep2'); output_prefix <- 'ATAC_MEF_Dox_D1'
#groups <- c('MEF_Dox_D2_rep1', 'MEF_Dox_D2_rep2'); output_prefix <- 'ATAC_MEF_Dox_D2'
#groups <- c('MEF_Dox_D7_Flk1pos_rep1', 'MEF_Dox_D7_Flk1pos_rep2'); output_prefix <- 'ATAC_MEF_Dox_D7_Flk1pos'
#groups <- c('EB_Dox_D25_rep1', 'EB_Dox_D25_rep2'); output_prefix <- 'ATAC_EB_Dox_D25'
groups <- c('EB_Dox_D25_rep1', 'EB_Dox_D25_rep2'); output_prefix <- 'ATAC_EB_Dox_D25'
input.bam.files <- d[d[, 'group'] %in% groups,'bam.file']
output.bam.file <- sprintf('analysis/etv2_pioneer/results/%s.bam', output_prefix)
source('sequencing.r'); merge.bam.files(input.bam.files, output.bam.file)
source('chipseq.r'); nucleoatac(output.bam.file, peak_file2, genome = 'mm10', mc.cores = 8)


# --------------------------------------------------------------------------
source('chipseq.r')
peak_file <- sprintf('analysis/etv2_pioneer/results/%s_summits.bed', 'MEF_Dox_D1_Etv2')
gr <- read.table(peak_file, header = FALSE, sep = '\t')
gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]), peak_id = gr[, 4], log10pvalue = -gr[, 5])
gr <- gr[values(gr)[['log10pvalue']] < -10]
gr <- resize(gr, width = 2000, fix = 'center')	# extending centered at the summit
gr <- add.seqinfo(gr, genome = 'mm10')
bw_files <- c(
	NoDox = sprintf('analysis/etv2_pioneer/results/ATAC_MEF_NoDox.nucleoatac.occ.bw'),
	Dox_D1 = sprintf('analysis/etv2_pioneer/results/ATAC_MEF_Dox_D1.nucleoatac.occ.bw'),
	Dox_D2 = sprintf('analysis/etv2_pioneer/results/ATAC_MEF_Dox_D2.nucleoatac.occ.bw'),
	Dox_D7 = sprintf('analysis/etv2_pioneer/results/ATAC_MEF_Dox_D7.nucleoatac.occ.bw')
)
cvgs <- lapply(bw_files, function(bw_file) import(bw_file , which = trim(reduce(gr)), as = 'RleList'))

Z <- lapply(cvgs, function(cvg) as(as(cvg[gr], 'RleViews'), 'matrix'))
Y <- do.call('rbind', lapply(Z, colMeans))

plot(Y[1, ], ylim = range(Y), type = 'l', lwd = 3, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
lines(Y[2, ], col = 'red', lwd = 3)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
pmt <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene)
is_promoter <- 1:length(gr) %in% as.matrix(findOverlaps(gr, pmt))[, 1]

H3K27ac_bw_file <- 'analysis/etv2_pioneer/results/MEF_NoDox_H3K27ac_treat_pileup.bw'
cvg_H3K27ac <- import(H3K27ac_bw_file, which = trim(reduce(gr)), as = 'RleList')
Z_H3K27ac <- as(as(cvg_H3K27ac[gr], 'RleViews'), 'matrix')
#p <- rowMeans(Z_H3K27ac[, 800:1200])	# the pileup density at the center [-100, +100] region
p <- rowMeans(Z_H3K27ac)	# the pileup density at the center [-100, +100] region
not_enhancer <- 1:length(gr) %in% order(p, decreasing = FALSE)[1:5000]
plot(colMeans(Z_H3K27ac[not_enhancer, ]))


#i <- !is_promoter & !is_enhancer
#i <- !is_promoter & not_enhancer
i <- not_enhancer
#plot(colMeans(Z[['NoDox']][i, ]), type = 'l', lwd = 3, ylim = c(0.5, 0.7), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
plot(colMeans(Z[['NoDox']][i, ]), type = 'l', lwd = 3, ylim = c(0.52, 0.58))
lines(colMeans(Z[['Dox_D1']][i, ]), lwd = 3, col = 'red')

i <- !not_enhancer
plot(colMeans(Z[['NoDox']][i, ]), type = 'l', lwd = 3, ylim = c(0.5, 0.7))
lines(colMeans(Z[['Dox_D1']][i, ]), lwd = 3, col = 'red')

cvg_brg1 <- import('analysis/etv2_pioneer/results/MEF_NoDox_Brg1_treat_pileup.bw', which = trim(reduce(gr)), as = 'RleList')
Z_Brg1_NoDox <- as(as(cvg_brg1[gr], 'RleViews'), 'matrix')
Z_Brg1_NoDox <- t(apply(Z_Brg1_NoDox, 1, function(z) (z - min(z)) / (max(z) - min(z))))
plot(colMeans(Z_Brg1_NoDox[i, ]), type = 'l', lwd = 3)


# --------------------------------------------------------------------------
# [2019-03-14] Generate two sets of Etv2 binding sites at D1
# 1. Etv2 peaks at H3K27ac
# 2. Etv2 peaks w/o H3K27ac
# [2019-05-07] Generate the epigenomic heatmap for more histone marks from Chronis et al.
# This is for preparing a slide for Ken's visit on 2019-05-08
# --------------------------------------------------------------------------
source('chipseq.r')
peak_file <- sprintf('analysis/etv2_pioneer/results/%s_summits.bed', 'MEF_Dox_D1_Etv2')
gr <- read.table(peak_file, header = FALSE, sep = '\t')
gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]), peak_id = gr[, 4], log10pvalue = -gr[, 5])
gr <- gr[values(gr)[['log10pvalue']] < -10]
gr <- resize(gr, width = 2000, fix = 'center')  # extending centered at the summit
gr <- add.seqinfo(gr, genome = 'mm10')

H3K27ac_bw_file <- 'analysis/etv2_pioneer/results/MEF_NoDox_H3K27ac_treat_pileup.bw'
cvg_H3K27ac <- import(H3K27ac_bw_file, which = trim(reduce(gr)), as = 'RleList')
Z_H3K27ac <- as(as(cvg_H3K27ac[gr], 'RleViews'), 'matrix')
p <- rowMeans(Z_H3K27ac)  

bw_file <- 'analysis/etv2_pioneer/data/Chronis/MNase_treat_pileup.bw'; high_color <- 'green'
bw_file <- 'analysis/etv2_pioneer/data/Chronis/H3K4me3_treat_pileup.bw'; high_color <- 'orange'
bw_file <- 'analysis/etv2_pioneer/data/Chronis/Brg1_treat_pileup.bw'; high_color <- 'red'
bw_file <- 'analysis/etv2_pioneer/data/Chronis/H3K36me3_treat_pileup.bw'; high_color <- 'brown'	# no clear signal
bw_file <- 'analysis/etv2_pioneer/data/Chronis/H3K9ac_treat_pileup.bw'; high_color <- 'hotpink'
bw_file <- 'analysis/etv2_pioneer/data/Chronis/H3K79me2_treat_pileup.bw'; high_color <- 'brown'
bw_file <- 'analysis/etv2_pioneer/data/Chronis/H3K27ac_treat_pileup.bw'; high_color <- 'blue'
bw_file <- 'analysis/etv2_pioneer/data/Chronis/H3K27me3_treat_pileup.bw'; high_color <- 'olivedrab1'
bw_file <- 'analysis/etv2_pioneer/data/Chronis/H3K4me2_treat_pileup.bw'; high_color <- 'mediumpurple'
bw_file <- 'analysis/etv2_pioneer/data/Chronis/H3K4me1_treat_pileup.bw'; high_color <- 'lightgoldenrod'
bw_file <- 'analysis/etv2_pioneer/data/Chronis/H3K9me3_treat_pileup.bw'; high_color <- 'lightblue'
bw_file <- 'analysis/etv2_pioneer/data/Chronis/Hdac1_treat_pileup.bw'; high_color <- 'orchid'
bw_file <- 'analysis/etv2_pioneer/data/Chronis/ATAC_treat_pileup.bw'; high_color <- 'skyblue'
bw_file <- 'analysis/etv2_pioneer/data/Chronis/H3.3_treat_pileup.bw'; high_color <- 'rosybrown1'
bw_file <- 'analysis/etv2_pioneer/data/Chronis/H3_treat_pileup.bw'; high_color <- 'yellow'
bw_file <- 'analysis/etv2_pioneer/data/Chronis/P300_treat_pileup.bw'; high_color <- 'gold'

cvg <- import(bw_file, which = trim(reduce(gr)), as = 'RleList')
Z <- as(as(cvg[gr], 'RleViews'), 'matrix')
image(t(Z[order(p), ]), breaks = quantile(Z, seq(0, 1, length.out = 101)), col = colorpanel(100, low = 'black', mid = 'gray', high = high_color), axes = FALSE, main = bw_file)


H3K27ac_D1_bw_file <- 'analysis/etv2_pioneer/results/MEF_Dox_D1_H3K27ac_treat_pileup.bw'
cvg_H3K27ac_D1 <- import(H3K27ac_D1_bw_file, which = trim(reduce(gr)), as = 'RleList')
Z_H3K27ac_D1 <- as(as(cvg_H3K27ac_D1[gr], 'RleViews'), 'matrix')

Etv2_bw_file <- 'analysis/etv2_pioneer/results/MEF_Dox_D1_Etv2_treat_pileup.bw'
cvg_Etv2 <- import(Etv2_bw_file, which = trim(reduce(gr)), as = 'RleList')
Z_Etv2 <- as(as(cvg_Etv2[gr], 'RleViews'), 'matrix')

Brg1_bw_file <- 'analysis/etv2_pioneer/results/MEF_NoDox_Brg1_treat_pileup.bw'
cvg_Brg1 <- import(Brg1_bw_file, which = trim(reduce(gr)), as = 'RleList')
Z_Brg1 <- as(as(cvg_Brg1[gr], 'RleViews'), 'matrix')

Brg1_D1_bw_file <- 'analysis/etv2_pioneer/results/MEF_Dox_D1_Brg1_treat_pileup.bw'
cvg_Brg1_D1 <- import(Brg1_D1_bw_file, which = trim(reduce(gr)), as = 'RleList')
Z_Brg1_D1 <- as(as(cvg_Brg1_D1[gr], 'RleViews'), 'matrix')


image(t(Z_H3K27ac[order(p), ]), breaks = quantile(Z_H3K27ac, seq(0, 1, length.out = 101)), col = colorpanel(100, low = 'black', mid = 'gray', high = 'blue'), axes = FALSE)
image(t(Z_H3K27ac_D1[order(p), ]), breaks = quantile(Z_H3K27ac_D1, seq(0, 1, length.out = 101)), col = colorpanel(100, low = 'black', mid = 'gray', high = 'blue'), axes = FALSE)
image(t(Z_Etv2[order(p), ]), breaks = quantile(Z_Etv2, seq(0, 1, length.out = 101)), col = colorpanel(100, low = 'black', mid = 'gray', high = 'white'), axes = FALSE)
image(t(Z_Brg1[order(p), ]), breaks = quantile(Z_Brg1, seq(0, 1, length.out = 101)), col = colorpanel(100, low = 'black', mid = 'gray', high = 'red'), axes = FALSE)
image(t(Z_Brg1_D1[order(p), ]), breaks = quantile(Z_Brg1_D1, seq(0, 1, length.out = 101)), col = colorpanel(100, low = 'black', mid = 'gray', high = 'red'), axes = FALSE)

#ATAC_NoDOX_bw_file <- 'analysis/etv2_pioneer/results/ATAC_MEF_NoDox.nucleoatac.occ.bw'
ATAC_NoDOX_bw_file <- 'analysis/etv2_pioneer/results/ATAC_MEF_NoDox.nucleoatac.nucleoatac_signal.smooth.bw'
#ATAC_NoDOX_bw_file <- 'analysis/etv2_pioneer/results/ATAC_MEF_NoDox_treat_pileup.bw'
cvg_ATAC_NoDOX <- import(ATAC_NoDOX_bw_file, which = trim(reduce(gr)), as = 'RleList')
Z_ATAC_NoDOX <- as(as(cvg_ATAC_NoDOX[gr], 'RleViews'), 'matrix')

#ATAC_D1_bw_file <- 'analysis/etv2_pioneer/results/ATAC_MEF_Dox_D1.nucleoatac.occ.bw'
#ATAC_D1_bw_file <- 'analysis/etv2_pioneer/results/ATAC_MEF_Dox_D1.nucleoatac.nucleoatac_signal.smooth.bw'
ATAC_D1_bw_file <- 'analysis/etv2_pioneer/results/ATAC_MEF_Dox_D1_treat_pileup.bw'
cvg_ATAC_D1 <- import(ATAC_D1_bw_file, which = trim(reduce(gr)), as = 'RleList')
Z_ATAC_D1 <- as(as(cvg_ATAC_D1[gr], 'RleViews'), 'matrix')

ATAC_D1_vs_NoDox_bw_file <- 'analysis/etv2_pioneer/results/ATAC_Dox_D1_vs_NoDox_treat_pileup.bw'
cvg_ATAC_D1_vs_NoDox_bw_file <- import(ATAC_D1_vs_NoDox_bw_file, which = trim(reduce(gr)), as = 'RleList')
Z_ATAC_D1_vs_NoDox <- as(as(cvg_ATAC_D1_vs_NoDox_bw_file[gr], 'RleViews'), 'matrix')


# --------------------------------------------------------------------------
# [2019-03-19] Compare and visualize the MEF MNase-seq signal of two groups of Etv2 peaks
# 1. Etv2 peaks with H3K27ac
# 2. Etv2 peaks w/o H3K27ac
# --------------------------------------------------------------------------
MNase_bw_file <- 'analysis/etv2_pioneer/data/MEF_Coverage_pileup.bw'
cvg_MNase <- import(MNase_bw_file, which = trim(reduce(gr)), as = 'RleList')
Z_MNase <- as(as(cvg_MNase[gr], 'RleViews'), 'matrix')

image(t(Z_ATAC_NoDOX[order(p), ]), breaks = quantile(Z_ATAC_NoDOX, seq(0, 1, length.out = 101)), col = colorpanel(100, low = 'black', mid = 'gray', high = 'green'), axes = FALSE)
image(t(Z_ATAC_D1[order(p), ]), breaks = quantile(Z_ATAC_D1, seq(0, 1, length.out = 101)), col = colorpanel(100, low = 'black', mid = 'gray', high = 'green'), axes = FALSE)

ZZ <- rbind(
	colMeans(Z_MNase[order(p, decreasing = FALSE)[1:10000], 800:1200]),
	colMeans(Z_MNase[order(p, decreasing = TRUE)[1:10000], 800:1200])
)
ZZ <- t(apply(ZZ, 1, function(zz) (zz - min(zz)) / (max(zz) - min(zz))))
plot(ZZ[1, ], xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', ylim = c(0, 1), lwd = 3, type = 'l', col = 'darkgreen')
lines(ZZ[2, ], col = 'darkgreen', lwd = 3, lty = 2)
abline(v = 200)




# ----------------------------------------------------------------------------
# [2019-03-14] Calling ATAC-seq peaks using MACS2 by comparing D1 and NoDox
# ----------------------------------------------------------------------------
dataset <- 'dataset=Etv2ATAC_version=20190228a'
source('sra.r'); d <- read.dataset(dataset, touch = FALSE)
treatment <- d[, 'group'] %in% c('MEF_Dox_D1_rep1', 'MEF_Dox_D1_rep2')
control <- d[, 'group'] %in% c('MEF_NoDox_rep1', 'MEF_NoDox_rep2')
base.name <- sprintf('analysis/etv2_pioneer/results/ATAC_Dox_D1_vs_NoDox')
pileup_file <- sprintf('%s_treat_pileup.bw', base.name)
s3_pileup_file <- gsub('.+/(.+)', 's3://etv2_pioneer/\\1', pileup_file)
s3_pileup_public_file <- gsub('.+/(.+)', 'https://s3.msi.umn.edu/etv2_pioneer/\\1', pileup_file)

treatment_files <- d[treatment, 'bam.file']
control_files <- d[control, 'bam.file']
source('chipseq.r'); macs2.callpeak(treatment_files, base.name, control_files, format = 'BAMPE', genome = 'mm10', broad = FALSE, qvalue.cutoff = 0.05, fold.change = FALSE, update = TRUE, call.summits = TRUE, shift = -100, extsize = 200)



# ----------------------------------------------------------------------------
# [2019-03-14] Calling Brg1 peaks using MACS2 by comparing D1 and NoDox
# ----------------------------------------------------------------------------
dataset <- 'dataset=Etv2ChIPseq_version=20190307a'
source('sra.r'); d <- read.dataset(dataset, touch = FALSE)
treatment <- d[, 'group'] %in% c('MEF_Dox_D1_Brg1_rep1', 'MEF_Dox_D1_Brg1_rep2')
control <- d[, 'group'] %in% c('MEF_NoDox_Brg1_rep1', 'MEF_NoDox_Brg1_rep2')
base.name <- sprintf('analysis/etv2_pioneer/results/MEF_Dox_D1_vs_NoDox_Brg1')
pileup_file <- sprintf('%s_treat_pileup.bw', base.name)

treatment_files <- d[treatment, 'bam.file']
control_files <- d[control, 'bam.file']
source('chipseq.r'); macs2.callpeak(treatment_files, base.name, control_files, format = 'BAM', genome = 'mm10', broad = FALSE, qvalue.cutoff = 0.05, fold.change = TRUE, update = TRUE, call.summits = TRUE)


# ----------------------------------------------------------------------------
# [2019-03-14] Look at the dynamics of Etv2 peaks
# ----------------------------------------------------------------------------
se_file <- sprintf('analysis/etv2_pioneer/results/MEF_Etv2_peaks.rds')
gr <- readRDS(se_file)
gr <- resize(gr, width = 2000, fix = 'center')	# extending centered at the summit

H3K27ac_bw_file <- 'analysis/etv2_pioneer/results/MEF_NoDox_H3K27ac_treat_pileup.bw'
cvg_H3K27ac <- import(H3K27ac_bw_file, which = trim(reduce(gr)), as = 'RleList')
Z_H3K27ac <- as(as(cvg_H3K27ac[gr], 'RleViews'), 'matrix')
p <- rowMeans(Z_H3K27ac)  

Etv2_bw_files <- c(
	'analysis/etv2_pioneer/results/MEF_Dox_D1_Etv2_treat_pileup.bw',
	'analysis/etv2_pioneer/results/MEF_Dox_D2_Etv2_treat_pileup.bw',
	'analysis/etv2_pioneer/results/MEF_Dox_D7_Etv2_treat_pileup.bw'
)
Z_Etv2 <- lapply(Etv2_bw_files, function(bw_file){
	cvg <- import(bw_file, which = trim(reduce(gr)), as = 'RleList')
	as(as(cvg[gr], 'RleViews'), 'matrix')
})
lapply(Z_Etv2, function(Z){
	image(t(Z[order(p), ]), breaks = quantile(Z, seq(0, 1, length.out = 101)), col = colorpanel(100, low = 'black', mid = 'gray', high = 'white'), axes = FALSE)
})


H3K27ac_bw_files <- c(
	'analysis/etv2_pioneer/results/MEF_NoDox_H3K27ac_treat_pileup.bw',
	'analysis/etv2_pioneer/results/MEF_Dox_D1_H3K27ac_treat_pileup.bw',
	'analysis/etv2_pioneer/results/MEF_Dox_D2_H3K27ac_treat_pileup.bw',
	'analysis/etv2_pioneer/results/MEF_Dox_D7_H3K27ac_treat_pileup.bw'
)
Z_H3K27ac <- lapply(H3K27ac_bw_files, function(bw_file){
	cvg <- import(bw_file, which = trim(reduce(gr)), as = 'RleList')
	as(as(cvg[gr], 'RleViews'), 'matrix')
})

lapply(Z_H3K27ac, function(Z){
	image(t(Z[order(p), ]), breaks = quantile(Z, seq(0, 1, length.out = 101)), col = colorpanel(100, low = 'black', mid = 'gray', high = 'blue'), axes = FALSE)
})
plot(colMeans(Z_H3K27ac[[1]]), col = 'blue', lwd = 3)


Brg1_bw_files <- c(
	'analysis/etv2_pioneer/results/MEF_NoDox_Brg1_treat_pileup.bw',
	'analysis/etv2_pioneer/results/MEF_Dox_D1_Brg1_treat_pileup.bw',
	'analysis/etv2_pioneer/results/MEF_Dox_D2_Brg1_treat_pileup.bw',
	'analysis/etv2_pioneer/results/MEF_Dox_D7_Brg1_treat_pileup.bw'
)
Z_Brg1 <- lapply(Brg1_bw_files, function(bw_file){
	cvg <- import(bw_file, which = trim(reduce(gr)), as = 'RleList')
	as(as(cvg[gr], 'RleViews'), 'matrix')
})

lapply(Z_Brg1, function(Z){
	image(t(Z[order(p), ]), breaks = quantile(Z, seq(0, 1, length.out = 101)), col = colorpanel(100, low = 'black', mid = 'gray', high = 'red'), axes = FALSE)
})



ATAC_bw_files <- c(
	'analysis/etv2_pioneer/results/ATAC_MEF_NoDox_treat_pileup.bw',
	'analysis/etv2_pioneer/results/ATAC_MEF_Dox_D1_treat_pileup.bw',
	'analysis/etv2_pioneer/results/ATAC_MEF_Dox_D2_treat_pileup.bw',
	'analysis/etv2_pioneer/results/ATAC_MEF_Dox_D7_treat_pileup.bw'
)
Z_ATAC <- lapply(ATAC_bw_files, function(bw_file){
	cvg <- import(bw_file, which = trim(reduce(gr)), as = 'RleList')
	as(as(cvg[gr], 'RleViews'), 'matrix')
})

lapply(Z_ATAC, function(Z){
	image(t(Z[order(p), ]), breaks = quantile(Z, seq(0, 1, length.out = 101)), col = colorpanel(100, low = 'black', mid = 'gray', high = 'green'), axes = FALSE)
})



# ----------------------------------------------------------------------------
# [2019-03-14] Difference of Brg1 after Etv2 induction
# [2019-04-18] Visualize anything's (e.g. Brg1) intensity surrounding Etv2 ChIP-seq peaks
# ----------------------------------------------------------------------------
source('chipseq.r')
peak_file <- sprintf('analysis/etv2_pioneer/results/%s_summits.bed', 'MEF_Dox_D1_Etv2')
gr <- read.table(peak_file, header = FALSE, sep = '\t')
gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]), peak_id = gr[, 4], log10pvalue = -gr[, 5])
gr <- gr[values(gr)[['log10pvalue']] < -10]
gr <- resize(gr, width = 2000, fix = 'center')  # extending centered at the summit
gr <- add.seqinfo(gr, genome = 'mm10')

# The H3K27ac signal in MEF; The Etv2 peaks will be sorted based on H3K27ac
H3K27ac_bw_file <- 'analysis/etv2_pioneer/results/MEF_NoDox_H3K27ac_treat_pileup.bw'
cvg_H3K27ac <- import(H3K27ac_bw_file, which = trim(reduce(gr)), as = 'RleList')
Z_H3K27ac <- as(as(cvg_H3K27ac[gr], 'RleViews'), 'matrix')
p <- rowMeans(Z_H3K27ac)  

#Brg1_diff_bw_file <- 'analysis/etv2_pioneer/results/MEF_Dox_D1_vs_NoDox_Brg1_treat_pileup.bw'
#bw_file <- 'analysis/etv2_pioneer/results/MEF_Brg1_GSE90893_treat_pileup.bw'
bw_file <- 'analysis/etv2_pioneer/results/MEF_Brg1_GSE71507_treat_pileup.bw'
cvg <- import(bw_file, which = trim(reduce(gr)), as = 'RleList')
Z <- as(as(cvg[gr], 'RleViews'), 'matrix')
image(t(Z[order(p), ]), breaks = quantile(Z, seq(0, 1, length.out = 101)), col = colorpanel(100, low = 'black', mid = 'gray', high = 'red'), axes = FALSE)


# ----------------------------------------------------------------------------
# [2019-05-26] Generate bed files for H3K27ac low/high Etv2 peaks at D1 MEF
# Still use Homer instead of ChIPpeakAnno, which generate erros when using TxDB as the annotation
# database; Homer is more stable. 
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); etv2 <- read_Etv2_peaks_in_D1_MEF(exclude_promoter = FALSE, exclude_exons = FALSE, MEF_H3K27ac = TRUE)
peak_file <- sprintf('%s/macs2/Etv2_high_H3K27ac.bed', PROJECT_DIR)
write.table(as.data.frame(etv2)[, 1:3], peak_file , sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
source('chipseq.r'); homer(peak_file, 'mm10', mc.cores = 4)
output_dir <- sprintf('%s_homer', peak_file)
source('s3.r'); s3.backup(output_dir, 's3://etv2_pioneer/homer/', make_public = TRUE)


source('analysis/etv2_pioneer/helper.r'); etv2 <- read_Etv2_peaks_in_D1_MEF(exclude_promoter = FALSE, exclude_exons = FALSE, MEF_H3K27ac = FALSE)
peak_file <- sprintf('%s/macs2/Etv2_low_H3K27ac.bed', PROJECT_DIR)
write.table(as.data.frame(etv2)[, 1:3], peak_file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
source('chipseq.r'); homer(peak_file, 'mm10', mc.cores = 4)
output_dir <- sprintf('%s_homer', peak_file)
source('s3.r'); s3.backup(output_dir, 's3://etv2_pioneer/homer/', make_public = TRUE)


# ----------------------------------------------------------------------------
# [2019-05-26] Compare the chromosonal regions of H3K27ac high/low Etv2 peaks
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); high <- read_Etv2_peaks_in_D1_MEF(exclude_promoter = FALSE, exclude_exons = FALSE, MEF_H3K27ac = TRUE)
aCR_high <- assignChromosomeRegion(high , nucleotideLevel = FALSE,  precedence = c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Exons", "Introns"),  TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)

source('analysis/etv2_pioneer/helper.r'); low <- read_Etv2_peaks_in_D1_MEF(exclude_promoter = FALSE, exclude_exons = FALSE, MEF_H3K27ac = FALSE)
aCR_low <- assignChromosomeRegion(low, nucleotideLevel = FALSE,  precedence = c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Exons", "Introns"),  TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)

x <- rbind(
	data.frame(group = 'low', region = names(aCR_low$percentage), percent = as.numeric(aCR_low$percentage)),
	data.frame(group = 'high', region = names(aCR_high$percentage), percent = as.numeric(aCR_high$percentage))
)
x <- transform(x, region = factor(region, c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Exons", "Introns")))
library(ggplot2)
ggplot(x, aes(x = region, y = percent, fill = group)) + geom_bar(stat = 'identity', position=position_dodge()) + coord_flip() + theme(axis.text=element_text(size=16, face = 'bold')) + scale_fill_manual(values=c('red','black'))


# ----------------------------------------------------------------------------
# [2019-03-18] A Venn diagram of the Etv2 peaks across three time points
# ----------------------------------------------------------------------------
library(SummarizedExperiment)
source('chipseq.r')
peaks_files <- sprintf('analysis/etv2_pioneer/results/%s_peaks.narrowPeak', c(
		'MEF_Dox_D1_Etv2', 
		'MEF_Dox_D2_Etv2',
		'MEF_Dox_D7_Etv2'
))
grl <- lapply(peaks_files, function(file){
	gr <- read.table(file, header = FALSE, sep = '\t')
	gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]), peak_id = gr[, 4], log10pvalue = -gr[, 5])
	gr <- gr[values(gr)[['log10pvalue']] < -10]
	gr <- add.seqinfo(gr, genome = 'mm10')
	gr
})

gr0 <- Reduce('union', grl)
x <- lapply(1:length(grl), function(i){
	as.matrix(findOverlaps(gr0, grl[[i]]))[, 1]
})
names(x) <- c('D1', 'D2', 'D7')
library(VennDiagram)
venn.plot <- venn.diagram(x, filename = 'analysis/etv2_pioneer/results/Venn_Etv2_ChIP-seq.tiff', euler.d = TRUE)


# ----------------------------------------------------------------------------
# [2019-03-19] Testing the hypothesis that if the H3K27ac low Etv2 peaks locate 
# near the endothelial genes
# ----------------------------------------------------------------------------
low_output_dir <- sprintf('%s/macs2/Etv2_low_H3K27ac.bed_homer', PROJECT_DIR)
low <- read.table(sprintf('%s/go/biological_process.txt', low_output_dir), header = TRUE, sep = '\t', comment.char = '', quote = '', fill = TRUE)

high_output_dir <- sprintf('%s/macs2/Etv2_high_H3K27ac.bed_homer', PROJECT_DIR)
high <- read.table(sprintf('%s/go/biological_process.txt', high_output_dir), header = TRUE, sep = '\t', comment.char = '', quote = '', fill = TRUE)

y <- merge(data.frame(term = low[, 2], go = low[, 1], p_low = low[, 3]), data.frame(go = high[, 1], p_high = high[, 3]), by.x = 'go', by.y = 'go')
y <- transform(y, score = -log((p_low + 1e-10) / (p_high + 1e-10)))

y_low <- y[y$p_low < 1e-3, ]
y_low <- y_low[order(y_low[, 'score'], decreasing = TRUE), ]
low_terms <- c('signaling', 'cell communication', 'ion transport', 'chemotaxis', 'blood circulation', 'mesenchymal cell development', 'stem cell differentiation', 'neuron projection guidance')
y_low <- y[y$term %in% low_terms, ]
y_low <- y_low[order(y_low$p_low), ]
y_low <- transform(y_low, term = factor(term, rev(y_low$term)))
library(ggplot2)
ggplot(data = y_low, aes(x = term, y = -log10(p_low))) + geom_bar(stat = 'identity') + coord_flip() + theme(axis.text=element_text(size=16, face = 'bold'))

y_high <- y[y$p_high < 1e-3, ]
y_high <- y_high[order(y_high[, 'score'], decreasing = FALSE), ]
high_terms <- c('translation', 'protein acetylation', 'nucleic acid transport', 'histone modification', 'regulation of RNA splicing', 'RNA localization', 'RNA catabolic process', 'protein folding')
y_high <- y[y$term %in% high_terms, ]
y_high <- y_high[order(y_high$p_high), ]
y_high <- transform(y_high , term = factor(term, rev(y_high$term)))
library(ggplot2)
ggplot(data = y_high, aes(x = term, y = -log10(p_high))) + geom_bar(stat = 'identity') + coord_flip() + theme(axis.text=element_text(size=16, face = 'bold'))



#i <- 'cell migration'
#i <- 'blood vessel development'
#i <- 'cell cycle'
#i <- 'chromatin organization'
i <- 'cell surface receptor signaling pathway'

low[which(low[, 2] == 'blood vessel development'), ]


# ----------------------------------------------------------------------------
# [2019-03-19]  Check if Etv2 binding is dependend on Brg1
# Results: The real Etv2 peaks have strong Brg1 signals, while the predicted Etv2 peaks
# have no Brg1 signals
# ----------------------------------------------------------------------------
source('chipseq.r')
peak_file <- sprintf('analysis/etv2_pioneer/results/%s_summits.bed', 'MEF_Dox_D1_Etv2')
gr <- read.table(peak_file, header = FALSE, sep = '\t')
gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]), peak_id = gr[, 4], log10pvalue = -gr[, 5])
gr <- resize(gr, width = 2000, fix = 'center')  # extending centered at the summit
real <- add.seqinfo(gr, genome = 'mm10')

library(chromVARmotifs) # https://github.com/GreenleafLab/chromVARmotifs
library(chromVAR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BiocParallel)
library(motifmatchr)
register(MulticoreParam(4)) # Use 8 cores
gr <- GRanges(seqnames = seqnames(BSgenome.Mmusculus.UCSC.mm10), range = IRanges(1, seqlengths(BSgenome.Mmusculus.UCSC.mm10)))
source('aux.r'); gr <- add.seqinfo(gr, genome = 'mm10')
se <- SummarizedExperiment(rowRanges = gr)
motif.set <- 'homer_pwms'
data(list = motif.set, package = 'chromVARmotifs')
motif_ix <- matchMotifs(get(motif.set)[82], se, genome = 'mm10', out = 'positions')
predicted <- motif_ix[[1]]
predicted <- etv2[seqnames(predicted) %in% sprintf('chr%s', c(1:19, 'X', 'Y'))]

mm <- as.matrix(findOverlaps(predicted, real))
neg <- predicted[!1:length(predicted) %in% mm[, 1]]
set.seed(1); neg <- neg[sample(1:length(neg), 50000)]
neg <- resize(neg, width = 2000, fix = 'center')


Brg1_bw_file <- 'analysis/etv2_pioneer/results/MEF_NoDox_Brg1_treat_pileup.bw'
cvg_Brg1 <- import(Brg1_bw_file, which = trim(reduce(c(real, neg))), as = 'RleList')

Z_neg <- as(as(cvg_Brg1[neg], 'RleViews'), 'matrix')
Z_real <- as(as(cvg_Brg1[real], 'RleViews'), 'matrix')

z_real <- colMeans(Z_real); z_real <- z_real - z_real[1]
z_neg <- colMeans(Z_neg); z_neg <- z_neg - z_neg[1]

plot(z_real, type = 'l', col = 'red', lwd = 3, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', ylim = range(z_real, z_neg))
lines(z_neg, col = 'gold', lwd = 3)


# ----------------------------------------------------------------------------
# [2019-03-19] Testing the hypothesis that Brg1 helps Etv2 bind to the non-enhancer
# region
# ----------------------------------------------------------------------------
low_peak_file <- 'analysis/etv2_pioneer/results/Etv2_low_H3K27ac.bed'
gr <- read.table(low_peak_file, header = FALSE, sep = '\t')
gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]))
gr <- resize(gr, width = 2000, fix = 'center')  # extending centered at the summit
low <- add.seqinfo(gr, genome = 'mm10')

high_peak_file <- 'analysis/etv2_pioneer/results/Etv2_high_H3K27ac.bed'
gr <- read.table(high_peak_file, header = FALSE, sep = '\t')
gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]))
gr <- resize(gr, width = 2000, fix = 'center')  # extending centered at the summit
high <- add.seqinfo(gr, genome = 'mm10')

bw_file <- sprintf('analysis/etv2_pioneer/results/MEF_Dox_D1_vs_NoDox_Brg1_FE.bw')
cvg <- import(bw_file, which = trim(reduce(c(high, low))), as = 'RleList')
Z_low <- as(as(cvg[low], 'RleViews'), 'matrix')
Z_high <- as(as(cvg[high], 'RleViews'), 'matrix')

bw_file <- sprintf('analysis/etv2_pioneer/results/ATAC_Dox_D1_vs_NoDox_FE.bw')
cvg <- import(bw_file, which = trim(reduce(c(high, low))), as = 'RleList')
Z_low <- as(as(cvg[low], 'RleViews'), 'matrix')
Z_high <- as(as(cvg[high], 'RleViews'), 'matrix')

plot(log2(colMeans(Z_low)), type = 'l', lwd = 3, col = 'red', xaxt = 'n', xlab = '', ylab = '', ylim = range(log2(colMeans(Z_low)), log2(colMeans(Z_high))))
lines(log2(colMeans(Z_high)), type = 'l', lwd = 3, col = 'blue')
abline(h = 0, lwd = 2, lty = 2)


# ----------------------------------------------------------------------------
# [2019-04-03] Investigate the ChIP-seq of OSKM on inducing OSKM in mouse MEF
# The ChIP-seq peaks are downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90893
# ----------------------------------------------------------------------------
source('aux.r')
library(rtracklayer)	# for liftover actions

ch <- import.chain(sprintf('%s/ucsc/mm9/mm9ToMm10.over.chain', Sys.getenv('SHARED')))	# have to be uncompressed file

peak_files <- c(
	Oct4 = sprintf('analysis/etv2_pioneer/data/GSE90893/Oct4_48hours.txt'),	
	Sox2 = sprintf('analysis/etv2_pioneer/data/GSE90893/Sox2_48hours.txt'),	
	Klf4 = sprintf('analysis/etv2_pioneer/data/GSE90893/Klf4_48hours.txt'),	
	Myc = sprintf('analysis/etv2_pioneer/data/GSE90893/Myc_48hours.txt'),	
	Fra = sprintf('analysis/etv2_pioneer/data/GSE90893/Fra_48hours.txt'),
	Cebpa = sprintf('analysis/etv2_pioneer/data/GSE90893/Cebpa_48hours.txt'),
	Cebpb = sprintf('analysis/etv2_pioneer/data/GSE90893/Cebpb_48hours.txt'),
	Runx1 = sprintf('analysis/etv2_pioneer/data/GSE90893/Runx1_48hours.txt')
#	Brg1 = sprintf('analysis/etv2_pioneer/data/GSE90893/Brg1_48hours.txt')
)

grs <- lapply(1:length(peak_files), function(i){
	cat(sprintf('processing %s\n', peak_files[i]))
	x <- read.table(peak_files[i], header = FALSE, sep = '\t')[, 1]	# the coordinates are in mm9
	chrom <- gsub('^(.+):(.+)-(.+)$', '\\1', x)
	from <- as.numeric(gsub('^(.+):(.+)-(.+)$', '\\2', x))
	to <- as.numeric(gsub('^(.+):(.+)-(.+)$', '\\3', x))
	x <- GRanges(seqnames = chrom, range = IRanges(from, to))
	x <- add.seqinfo(x, genome = 'mm9')
	x <- trim(x)
	seqlevelsStyle(x) <- 'UCSC'  # necessary
	x <- unlist(liftOver(x, ch))
	x <- add.seqinfo(x, genome = 'mm10')
	peaks <- trim(resize(x, width = 2000, fix = 'center'))# 1000 bp regions surrounding the mononucleosomes
	peaks
})
names(grs) <- names(peak_files)

peak_file <- sprintf('analysis/etv2_pioneer/results/%s_summits.bed', 'MEF_Dox_D1_Etv2')
gr <- read.table(peak_file, header = FALSE, sep = '\t')
gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]), peak_id = gr[, 4], log10pvalue = -gr[, 5])
gr <- gr[values(gr)[['log10pvalue']] < -20]
gr <- resize(gr, width = 2000, fix = 'center')  # extending centered at the summit
grs[['Etv2']] <- gr

peak_file <- 'analysis/etv2_pioneer/data/GSE43916/GSE43916_Ascl1_MEF_48hrs_summits.bed'
gr <- read.table(peak_file, header = FALSE, sep = '\t')
gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]))
gr <- add.seqinfo(gr, genome = 'mm9')
gr <- trim(gr)
seqlevelsStyle(gr) <- 'UCSC'  # necessary
gr <- unlist(liftOver(gr, ch))
gr <- add.seqinfo(gr, genome = 'mm10')
gr <- trim(resize(gr, width = 2000, fix = 'center'))# 1000 bp regions surrounding the mononucleosomes
grs[['Ascl1']] <- gr

#Brg1_bw_file <- 'analysis/etv2_pioneer/results/MEF_NoDox_Brg1_treat_pileup.bw'	# our Brg1 ChIP-seq in MEF
#Brg1_bw_file <- 'analysis/etv2_pioneer/results/MEF_Brg1_GSE90893_treat_pileup.bw'		# Brg1 ChIP-seq from GSE90893
Brg1_bw_file <- 'analysis/etv2_pioneer/results/MEF_Brg1_GSE71507_treat_pileup.bw'		# Brg1 ChIP-seq from GSE71507
Y <- do.call('rbind', lapply(1:length(grs), function(i){
	cat(sprintf('processing %s\n', names(grs)[i]))
	gr <- grs[[i]]
	cvg_Brg1 <- import(Brg1_bw_file, which = trim(reduce(gr)), as = 'RleList')
	y <- as(as(cvg_Brg1[gr], 'RleViews'), 'matrix')
	y_mean <- colMeans(y); y_mean <- y_mean - y_mean[1]
	y_mean
}))
rownames(Y) <- names(grs)
Y <- Y[c("Fra", "Cebpa", "Cebpb", "Oct4","Sox2", "Klf4", "Myc", 'Ascl1', "Etv2"), ]
image(t(Y), col = colorpanel(100, low = 'black', mid = 'red', high = 'yellow'), axes = FALSE)

plot(y_mean)


# ----------------------------------------------------------------------------
# [2019-04-07] Calling Brg1 peaks in MEF using MACS2 
# This is the Brg1 ChIP-seq data downlaoded form https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90893
# ----------------------------------------------------------------------------
dataset <- 'dataset=Chronis_version=20190405a'
source('sra.r'); d <- read.dataset(dataset, touch = TRUE)
d <- transform(d, sra.run.dir = sra.run.dir(run), bam.file = sprintf('%s/%s.dedup.bam', sra.run.result.dir(run), run))
treatment <- d[, 'group'] %in% c('MEF_Brg1')
control <- d[, 'group'] %in% c('MEF_input')
base.name <- sprintf('analysis/etv2_pioneer/results/MEF_Brg1_GSE90893')
pileup_file <- sprintf('%s_treat_pileup.bw', base.name)

treatment_files <- d[treatment, 'bam.file']
control_files <- d[control, 'bam.file']
source('chipseq.r'); macs2.callpeak(treatment_files, base.name, control_files, format = 'BAM', genome = 'mm10', broad = FALSE, qvalue.cutoff = 0.05, fold.change = FALSE , update = TRUE, call.summits = TRUE)


# ----------------------------------------------------------------------------
# [2019-04-08] Calling Brg1 peaks in MEF using MACS2 
# This is the Brg1 ChIP-seq data downlaoded form https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71507
# ----------------------------------------------------------------------------
dataset <- 'dataset=Alver_version=20190407a'
source('sra.r'); d <- read.dataset(dataset, touch = TRUE)
d <- transform(d, sra.run.dir = sra.run.dir(run), bam.file = sprintf('%s/%s.dedup.bam', sra.run.result.dir(run), run))
treatment <- d[, 'group'] %in% c('MEF_Brg1')
control <- d[, 'group'] %in% c('MEF_input')
base.name <- sprintf('analysis/etv2_pioneer/results/MEF_Brg1_GSE71507')
pileup_file <- sprintf('%s_treat_pileup.bw', base.name)

treatment_files <- d[treatment, 'bam.file']
control_files <- d[control, 'bam.file']
source('chipseq.r'); macs2.callpeak(treatment_files, base.name, control_files, format = 'BAM', genome = 'mm10', broad = FALSE, qvalue.cutoff = 0.05, fold.change = FALSE , update = TRUE, call.summits = TRUE)


# ----------------------------------------------------------------------------
# [2019-04-12] Uploading two public Brg1 ChIP-seq tracks
# ----------------------------------------------------------------------------
source('s3.r');s3.backup('analysis/etv2_pioneer/results/MEF_Brg1_GSE90893_treat_pileup.bw', 's3://etv2_pioneer/data/MEF_Brg1_GSE90893_treat_pileup.bw', make_public = TRUE)
source('s3.r');s3.backup('analysis/etv2_pioneer/results/MEF_Brg1_GSE71507_treat_pileup.bw', 's3://etv2_pioneer/data/MEF_Brg1_GSE71507_treat_pileup.bw', make_public = TRUE)


# ----------------------------------------------------------------------------
# [2019-04-12] Comparing two Brg1 ChIP-seq peaks and see the consistency between them
# If they are largely consistent, this will reduce the difficulities of re-doing the Brg1 ChIP-seq
# ----------------------------------------------------------------------------
library(GenomicRanges); 
source('chipseq.r')
peak1_file <- 'analysis/etv2_pioneer/results/MEF_Brg1_GSE71507_summits.bed'
peak2_file <- 'analysis/etv2_pioneer/results/MEF_Brg1_GSE90893_summits.bed'
gr1 <- read.table(peak1_file, header = FALSE, sep = '\t')
gr2 <- read.table(peak2_file, header = FALSE, sep = '\t')
gr1 <- gr1[gr1[, 5] > 10, ]
gr2 <- gr2[gr2[, 5] > 10, ]
gr1 <- GRanges(seqnames = gr1[, 1], range = IRanges(gr1[, 2], gr1[, 3]))
gr2 <- GRanges(seqnames = gr2[, 1], range = IRanges(gr2[, 2], gr2[, 3]))
gr1 <- resize(gr1, width = 500, fix = 'center')  # extending centered at the summit
gr2 <- resize(gr2, width = 500, fix = 'center')  # extending centered at the summit
gr1 <- add.seqinfo(gr1, genome = 'mm10')
gr2 <- add.seqinfo(gr2, genome = 'mm10')
mm <- as.matrix(findOverlaps(gr1, gr2))
table(1:length(gr1) %in% mm[, 1])
table(1:length(gr2) %in% mm[, 2])


bw_files <- c(
	'analysis/etv2_pioneer/results/MEF_NoDox_Brg1_treat_pileup.bw',	# our Brg1 ChIP-seq in MEF
	'analysis/etv2_pioneer/results/MEF_Brg1_GSE90893_treat_pileup.bw',		# Brg1 ChIP-seq from GSE90893
	'analysis/etv2_pioneer/results/MEF_Brg1_GSE71507_treat_pileup.bw'		# Brg1 ChIP-seq from GSE71507
)
Y <- do.call('cbind', lapply(bw_files, function(bw_file){
	cat(sprintf('processing %s\n', bw_file))
	cvg <- import(bw_file, which = trim(reduce(gr2)), as = 'RleList')
	mean(cvg[gr2])
}))


# ----------------------------------------------------------------------------
# [2019-04-16] Preprocessing the Etv2 scRNA-seq data
# Need to run on lab queue, will fail on the k40
# ----------------------------------------------------------------------------
dataset <- 'dataset=Etv2scRNAseq_version=20190416a'
se <- readRDS(sprintf('analysis/datasets/%s.rds', dataset))
library(BiocParallel)
register(MulticoreParam(4))
source('sctools.r'); se <- scran.preprocess(se)
se2 <- as(se, 'SummarizedExperiment')
rowData(se2) <- rowData(se)
dataset2 <- 'dataset=Etv2scRNAseq_version=20190416b'
saveRDS(se2, file = sprintf('analysis/datasets/%s.rds', dataset2))


# ----------------------------------------------------------------------------
# [2019-04-17] scVI analysis of the scRNA-seq of Etv2 induction in mouse
# Run this on k40
# [2019-04-23] SCRAN gives less HVGs
# ----------------------------------------------------------------------------
library(SummarizedExperiment)
dataset <- 'dataset=Etv2scRNAseq_version=20190416b'
se <- readRDS(sprintf('analysis/datasets/%s.rds', dataset))
n <- rowData(se)$FDR < 0.05
latent <- 10; dataset3 <- sprintf('dataset=Etv2scRNAseq_latent=%d_version=20190416c', latent)
#latent <- 20; dataset3 <- sprintf('dataset=Etv2scRNAseq_latent=%d_version=20190416c', latent)
devtools::load_all('analysis/ias/packages/RscVI'); se2 <- VAE(se[n, ], n_latent = latent)
colData(se)$latent <- colData(se2)$latent
saveRDS(se, file = sprintf('analysis/datasets/%s.rds', dataset3))



# ----------------------------------------------------------------------------
# [2019-04-17] Analyze the data by VAE
# ----------------------------------------------------------------------------
library(Seurat)
library(SummarizedExperiment)
dataset <- 'dataset=Etv2scRNAseq_version=20190416a'
se <- readRDS(sprintf('analysis/datasets/%s.rds', dataset))
source('sctools.r'); set.seed(1); se <- seurat.preprocess(se)
n_latent <- 10; dataset3 <- sprintf('dataset=Etv2scRNAseq_latent=%d_version=20190416d', n_latent)  # the filtered genes/cells after Seurat
devtools::load_all('analysis/ias/packages/RscVI'); se2 <- VAE(se, n_latent = n_latent)
saveRDS(se2, file = sprintf('analysis/datasets/%s.rds', dataset3))


# ----------------------------------------------------------------------------
# [2019-04-23] TSNE of the Etv2 scRNA-seq on MEF
# TSNE runs on the latent space obtained by scVI
# ----------------------------------------------------------------------------
library(SummarizedExperiment)
library(SingleCellExperiment)
latent <- 10; dataset <- sprintf('dataset=Etv2scRNAseq_latent=%d_version=20190416c', latent)
se <- readRDS(sprintf('analysis/datasets/%s.rds', dataset))

library(Rtsne); set.seed(1); y_tsne <- Rtsne(colData(se)$latent)$Y
group2bg <- c(
	'MEF_Dox_D1' = 'black', 
	'MEF_NoDox' = 'blue', 
	'MEF_Dox_D2' = 'purple', 
	'MEF_Dox_D7a' = 'red', 
	'MEF_Dox_D7b' = 'pink'
)
bg <- group2bg[colData(se)$group]
plot(y_tsne, cex = 0.5, pch = 21, bg  = bg, col = bg, xaxt = 'n', yaxt = 'n', main = 'cell2vec', xlab = '', ylab = '')


# ----------------------------------------------------------------------------
# [2019-04-25] Infer the pseudotime by slingshot
# clustering of cells on the TSNE space using Mclust
# It appears that the results of Mclust depend on the seed
# Mclust is used by slingshot and MST (see the Nat Biotechnol paper on comparing 
# different TI tools)
# [2019-04-05] examine the relationship cluster # and time points
# The cluster that has the most cells from NoDox will be used as the start
# [2019-04-27] Adding the cluster # to the t-SNE plot
# [2019-04-28] The Mclust results are not stable regarding the seeds; Use Louvain clustering instead
# ----------------------------------------------------------------------------
G_min <- 15; G_max <- 15 
library(mclust); set.seed(5); mc <- Mclust(y_tsne, G = G_min:G_max)
colData(se)$GMM <- mc$classification
library(RColorBrewer); plot(y_tsne, col = colorRampPalette(brewer.pal(11,'Spectral'))(G_max)[mc$classification], pch = 16, asp = 1, xaxt = 'n', yaxt = 'n', main = 'cluster', xlab = '', ylab = '')
set.seed(1); y_centers <- do.call('rbind', lapply(1:G_max, function(i) y_tsne[sample(which(mc$classification == i), 1), ]))
text(y_centers[, 1], y_centers[, 2], 1:G_max, cex = 3)
table(colData(se)$group, mc$classification)


# ----------------------------------------------------------------------------
# [2019-04-28] KNN graph and Louvain clustering
# ----------------------------------------------------------------------------
library(igraph)
library(FNN)
k <- 200
knn <- get.knn(y_tsne, k = k)
knn <- data.frame(from = rep(1:nrow(knn$nn.index), k), to = as.vector(knn$nn.index), weight = 1/(1 + as.vector(knn$nn.dist)))
g <- graph_from_data_frame(knn, directed = FALSE)
g <- simplify(g)
lc <- cluster_louvain(g)
clust <- as.numeric(as.factor(membership(lc)))
G_max <- max(unique(clust))
library(RColorBrewer); plot(y_tsne, col = colorRampPalette(brewer.pal(11,'Spectral'))(G_max)[clust], pch = 16, asp = 1, xaxt = 'n', yaxt = 'n', main = 'cluster', xlab = '', ylab = '')
set.seed(1); y_centers <- do.call('rbind', lapply(1:G_max, function(i) y_tsne[sample(which(clust == i), 1), ]))
text(y_centers[, 1], y_centers[, 2], 1:G_max, cex = 3)
table(colData(se)$group, clust)


# ----------------------------------------------------------------------------
# [2019-04-25] Use slingshot to get the lineages (a MST on the Mclust clusters)
# ----------------------------------------------------------------------------
start.clus <- '7'
library(slingshot); set.seed(1); lin <- getLineages(y_tsne, clust, start.clus = start.clus)
plot(y_tsne, col = bg, asp = 1, pch = 16, xaxt = 'n', yaxt = 'n', main = 'Lineages', xlab = '', ylab = '')
lines(lin, lwd = 3, show.constraints = TRUE)

library(RColorBrewer); plot(y_tsne, col = colorRampPalette(brewer.pal(11,'Spectral'))(G_max)[clust], pch = 16, asp = 1, xaxt = 'n', yaxt = 'n', main = 'cluster', xlab = '', ylab = '')
lines(lin, lwd = 3, show.constraints = TRUE)
text(y_centers[, 1], y_centers[, 2], 1:G_max, cex = 3)


# ----------------------------------------------------------------------------
# [2019-04-25] Get the smooth curve for each lineage
# This step takes roughly ~2 hours; the intermediate results are saved and re-used later
# [2019-04-30] Setting extend to 'n' is critical for stablize the curve, otherwise some 
# curve may go back to the start points and form a loop
# ----------------------------------------------------------------------------
set.seed(1); crv <- getCurves(lin, extend = 'n')	# this step takes a long time
crv_file <- sprintf('analysis/etv2_pioneer/results/scRNA-seq_Etv2_MEF_slingshot_curve.rds')
saveRDS(crv, file = crv_file)


crv_file <- sprintf('analysis/etv2_pioneer/results/scRNA-seq_Etv2_MEF_slingshot_curve.rds')
crv <- readRDS(crv_file)
plot(y_tsne, col = bg, asp = 1, pch = 16, xaxt = 'n', yaxt = 'n', main = 'trajectory', xlab = '', ylab = '')
library(slingshot)
lines(crv, lwd = 3, show.constraints = TRUE)


# ----------------------------------------------------------------------------
# [2019-04-28] Visualize gene expression levels by density plot
# ----------------------------------------------------------------------------
library(roxygen2); library(devtools); devtools::document('packages/denviz')


X <- assays(se)$logcounts; rownames(X) <- rowData(se)$name
devtools::load_all('packages/denviz'); denviz(X, y_tsne, g = 'Etv2', grid_points = 20)


# ----------------------------------------------------------------------------
# [2019-04-30] Get the density plot for dynamic genes
# ----------------------------------------------------------------------------
library(BiocParallel)
register(MulticoreParam(4)) # Use 8 cores
library(futile.logger); flog.threshold(TRACE)
X <- assays(se)$logcounts; rownames(X) <- rowData(se)$name; 
n <- rowData(se)$FDR < 0.05
grid_points <- 200
devtools::load_all('packages/denviz'); Z <- denmap(X[n, ], y_tsne, grid_points = grid_points)
rownames(Z) <- rownames(X)[n]
Z_file <- sprintf('analysis/etv2_pioneer/results/gene_cluster_grid=%d.rds', grid_points)
saveRDS(Z, Z_file)


# ----------------------------------------------------------------------------
# [2019-04-30] Clustering the genes based on their density
# ----------------------------------------------------------------------------
library(igraph)
library(FNN)
library(irlba); V <- prcomp_irlba(t(Z), n = 20)$rotation
k <- 50
grid_points <- 200
knn <- get.knn(Z, k = k)
knn <- data.frame(from = rep(1:nrow(knn$nn.index), k), to = as.vector(knn$nn.index), weight = 1/(1 + as.vector(knn$nn.dist)))
g <- graph_from_data_frame(knn, directed = FALSE)
g <- simplify(g)
lc <- cluster_louvain(g)
clust <- as.numeric(as.factor(membership(lc)))
names(clust) <- rownames(Z)

clust2 <- rep(NA, nrow(se))
clust2[n] <- clust
rowData(se)$cluster <- clust2
se_file <- sprintf('analysis/etv2_pioneer/data/%s_grid=%d_knn=%d.rds', dataset, grid_points, k)
saveRDS(se, se_file)


# ----------------------------------------------------------------------------
# [2019-04-30] t-SNE view of genes
# ----------------------------------------------------------------------------
grid_points <- 200; Z_file <- sprintf('analysis/etv2_pioneer/results/gene_cluster_grid=%d.rds', grid_points)
Z <- readRDS(Z_file)

library(irlba); V <- prcomp_irlba(t(Z), n = 20)$rotation
library(Rtsne); set.seed(1); g_tsne <- Rtsne(V)$Y
library(igraph)
library(FNN)
k <- 200
knn <- get.knn(g_tsne, k = k)
knn <- data.frame(from = rep(1:nrow(knn$nn.index), k), to = as.vector(knn$nn.index), weight = 1/(1 + as.vector(knn$nn.dist)))
g <- graph_from_data_frame(knn, directed = FALSE)
g <- simplify(g)
lc <- cluster_louvain(g)
clust <- as.numeric(as.factor(membership(lc)))
names(clust) <- rownames(Z)


grid_points <- 200; k <- 50
se_file <- sprintf('analysis/etv2_pioneer/data/%s_grid=%d_knn=%d.rds', dataset, grid_points, k)
se <- readRDS(se_file)

library(RColorBrewer); plot(g_tsne, col = colorRampPalette(brewer.pal(11,'Spectral'))(max(clust))[clust], pch = 16, asp = 1,, xaxt = 'n', yaxt = 'n', main = 'genes', xlab = '', ylab = '')
y_centers <- do.call('rbind', lapply(1:max(clust), function(i) c(median(g_tsne[clust == i, 1]), median(g_tsne[clust == i, 2]))))
text(y_centers[, 1], y_centers[, 2], 1:max(clust), cex = 3)



par(mfrow = c(3, 4), mar = c(0.5, 0.5, 0.5, 0.5))
lapply(1:max(clust), function(k){
Y <- matrix(colMeans(Z[clust == k, ]),grid_points,grid_points)
breaks <- c(seq(0, 0.001, length.out = 100), 1)
image(Y, col = colorpanel(100, low = 'black', mid = 'green', high = 'yellow'), breaks = breaks, xaxt = 'n', yaxt = 'n')
})


# ----------------------------------------------------------------------------
# [2019-04-30] Functional annotation of genes in each cluster
# Note done yet
# ----------------------------------------------------------------------------
grid_points <- 200; k <- 100
se_file <- sprintf('analysis/etv2_pioneer/data/%s_grid=%d_knn=%d.rds', dataset, grid_points, k)
se <- readRDS(se_file)
clust <- rowData(se)$cluster
names(clust) <- rowData(se)$name
clust <- clust[!is.na(clust)]
#lapply(split(clust, list(clust)), function(x) factor(names)

homer_file <- sprintf('analysis/etv2_pioneer/results/MEF_Dox_D7_Etv2_peaks.narrowPeak_homer/annotatePeaks.tsv')
x <- read.table(homer_file, sep = '\t', header = TRUE, comment.char = '', quote = '')
colnames(x)[1] <- 'peak'
annot <- gsub('(.+) \\(.+\\)', '\\1', x[, 'Annotation'])
annot <- gsub('\\.\\d+$','',  annot)
x <- cbind(x, Annotation2 = annot)

#x <- x[x[, 'Annotation2'] %in% c('promoter-TSS'), ]
x <- x[x[, 'Annotation2'] %in% c('Intergenic'), ]; x <- x[abs(x[, 'Distance.to.TSS']) < 1e5, ]

y <- unique(merge(data.frame(symbol = x[, 'Gene.Name']), data.frame(symbol = names(clust), clust = clust), by.x = 'symbol', by.y = 'symbol'))
barplot(table(y[, 2]) / table(clust))







# ----------------------------------------------------------------------------
# [2019-05-17] Plot expreession along a specified trajectory
# ----------------------------------------------------------------------------
X <- assays(se)$logcounts; rownames(X) <- rowData(se)$name

g <- 'Meis1'
x <- X[g, ]
endothelial <- 2	# the endothelial lineage
# get the cells and their order along the lineage #1
P <- slingPseudotime(crv)
n <- order(P[, endothelial])	# ordering cells based on their positions in lineage #1
n <- n[!is.na(P[n, endothelial])]
df <- data.frame(x = 1:length(n), y = x[n])
plot(y ~ x, df, pch = 21, bg = bg[n], col = bg[n], xlab = '', ylab = '', main = g)
#m <- df[, 'y'] > 0
m <- 1:nrow(df)
ss <- smooth.spline(df[m, 'x'], df[m, 'y'], df = 10)
lines(ss, lty = 2, col = 'black', lwd = 3)



sc <- slingCurves(crv)[[endothelial]]

points(sc$s, lwd = 2, col = 'green')


#

# Find the significant changed genes along the lineage #1
gs <- which(Matrix::rowSums(assays(se)$counts > 0) > 100)
res <- do.call('rbind', mclapply(gs, function(g){
	cat(sprintf('%d\n', g))
	x <- assays(se)$logcounts[g, ]
	df <- data.frame(x = 1:length(n), y = x[n] > 0)
	coef(summary(lm(y ~., df)))[2, ]
}, mc.cores = 4))


low_at_MEF <- Matrix::rowSums(assays(se)$counts[, colData(se)$group == 'MEF_Dox_D1']) / sum(colData(se)$group == 'MEF_Dox_D1') < 0.05
abundant <- Matrix::rowSums(assays(se)$counts[, n] > 0) > 20
gs <- which(low_at_MEF & abundant)

low_at_Flk1 <- Matrix::rowSums(assays(se)$counts[, colData(se)$group == 'MEF_Dox_D7b']) / sum(colData(se)$group == 'MEF_Dox_D7b') < 0.05
abundant <- Matrix::rowSums(assays(se)$counts[, n] > 0) > 20
gs <- which(low_at_Flk1 & abundant)

cc <- c(cor(t(as.matrix(assays(se)$logcounts[gs, n])), 1:length(n)))
names(cc) <- rowData(se)$name[gs]

#heatmap of top genes enriched at Flk1+ cell populations
gs2 <- names(sort(cc, decreasing = TRUE)[1:50])
#gs2 <- names(sort(cc, decreasing = FALSE)[1:50])

X <- assays(se)$logcounts; rownames(X) <- rowData(se)$name
X <- as.matrix(X[gs2, n])
colnames(X) <- sprintf('C%d', n)
breaks <- c(0, quantile(X[X > 0], seq(0, 1, length.out = 100)))
library(gplots); col2 <- colorpanel(100, low = 'black', mid = 'green', high = 'yellow')
df <- data.frame(time = colData(se)$group[n])
rownames(df) <- colnames(X)
collist <- list(time = group2bg)
library(pheatmap); pheatmap(X, cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = FALSE, breaks = breaks, color = col2, cellheight = 15, cellwidth = 0.03, fontsize_row = 15, show_colnames = FALSE, border_color = 'black', annotation_names_row = FALSE, annotation_names_col = FALSE, annotation_col = df, annotation_colors = collist)



# down-regulated
#gs2 <- names(sort(cc, decreasing = FALSE)[1:200])
gs2 <- names(sort(cc, decreasing = TRUE)[1:200])
z <- c('Cell Cycle' = 2.050E-46, 'Signaling by Rho GTPases' = 7.225E-15, 'PLK1 signaling events' = 8.743E-14, 'Aurora B signaling' = 6.044E-13, 'G1/S Transition' = 5.223E-8)
z <- c('vasculature development' = 1.035E-5, 'extracellular matrix organization' = 6.110E-5	, 'cardiovascular system development' = 2.389E-4, 'angiogenesis' = 3.142E-4	, 'endothelial cell migration' = 6.571E-3)


# ----------------------------------------------------------------------------
# [2019-04-28] Gene expression density plot
# ----------------------------------------------------------------------------
g <- 'Atf6'
h <- Matrix::colSums(assays(se)$logcounts) / nrow(se)
x <- assays(se)$logcounts[rowData(se)$name == g, ]
y <- y_tsne
n <- sample.int(ncol(se), 10000, prob = x / h, replace = TRUE)
lims <- c(range(y[, 1]) * 1.05, range(y[, 2]) * 1.05)
library(MASS); kd <- kde2d(y[n, 1], y[n, 2], n = 200, lims = lims)
library(gplots); col <- colorpanel(100, low = 'black', mid = 'green', high = 'yellow')
breaks <- c(seq(0, 0.001, length.out = 100), 1)
image(kd,  col = col, breaks = breaks, xaxt = 'n', yaxt = 'n',  main = g)


# ----------------------------------------------------------------------------
# [2019-04-28] Pathway analysis of gene expression on a t-SNE plot
# ----------------------------------------------------------------------------
	library(gplots); 
	library(MASS); 
se <- se[!duplicated(rowData(se)$name)]	# there exists duplicated gene names
source('go.r'); G <- go.matrix(rowData(se)$name, filter = 'external_gene_name', domain = 'BP', genome = 'mm10', method = 'ncbi')
ng <- Matrix::colSums(G)
G <- G[, ng < 2000 & ng > 10]

i <- 'GO:0001525'
X <- assays(se)$logcounts; rownames(X) <- rowData(se)$name
X <- X[rowSums(X > 0) > 20, ]
lapply(1:ncol(G), function(i){
	y <- y_tsne
	x <- colSums(X[G[, i], ])
	x <- X['Emcn', ]
	x <- rep(1, ncol(X))
	r <- diff(range(y)) / 100
	y <- y + rnorm(prod(dim(y)), mean = 0, sd = 0)
	n <- sample(1:length(x), 10000, prob = x, replace = TRUE)
	kd <- kde2d(y[n, 1], y[n, 2], n = 200, lims = c(range(y[, 1]) * 1.05, range(y[, 2]) * 1.05))
	col <- colorpanel(100, low = 'black', mid = 'green', high = 'yellow')
	image(kd,  col = col, xaxt = 'n', yaxt = 'n',  main = i)
})


# ----------------------------------------------------------------------------
# [2019-04-17] MACS2 on various MEF epigenic data
# ----------------------------------------------------------------------------
dataset <- 'dataset=Chronis_version=20190424a'
d <- read.table(sprintf('analysis/datasets/%s.tsv', dataset), sep = '\t', header = TRUE)
source('sra.r'); d <- transform(d, bam.file = sra.run.alignment.bam.file(run))

# Call narrow peaks
treatments <- c('ATAC', 'H3.3', 'P300', 'H3', 'Brg1', 'Hdac1')
source('chipseq.r');  library(parallel); mclapply(1:length(treatments), function(i){
	treatment_files <- d[d[, 'seqtype'] == treatments[i], 'bam.file']
	control_files <- d[d[, 'seqtype'] == 'input', 'bam.file']
	base.name <- sprintf('analysis/etv2_pioneer/data/Chronis/%s', treatments[i])
	macs2.callpeak(treatment_files, base.name, control_files, format = 'BAM', genome = 'mm10', broad = FALSE, qvalue.cutoff = 0.05, fold.change = FALSE , update = TRUE, call.summits = TRUE)
}, mc.cores = 4)

# Call broad peaks for histone modification
treatments <- c('H3K27me3', 'H3K4me1', 'H3K4me2', 'H3K9me3', 'H3K79me2', 'H3K36me3', 'H3K9ac', 'H3K27ac', 'H3K4me3')
source('chipseq.r');  library(parallel); mclapply(1:length(treatments), function(i){
	treatment_files <- d[d[, 'seqtype'] == treatments[i], 'bam.file']
	control_files <- d[d[, 'seqtype'] == 'input', 'bam.file']
	base.name <- sprintf('analysis/etv2_pioneer/data/Chronis/%s', treatments[i])
	macs2.callpeak(treatment_files, base.name, control_files, format = 'BAM', genome = 'mm10', broad = TRUE, qvalue.cutoff = 0.05, fold.change = FALSE , update = TRUE, call.summits = FALSE)
}, mc.cores = 4)


# Call peaks for MNase-seq
treatment <- c('MNase')
treatment_files <- d[d[, 'seqtype'] == treatment, 'bam.file']
base.name <- sprintf('analysis/etv2_pioneer/data/Chronis/%s', treatment)
source('chipseq.r'); macs2.callpeak(treatment_files, base.name, format = 'BAM', genome = 'mm10', broad = FALSE, pvalue.cutoff = 0.05, fold.change = FALSE , update = TRUE, call.summits = TRUE)


# ----------------------------------------------------------------------------
# [2019-04-17] MACS2 on various MEF epigenic data
# ----------------------------------------------------------------------------
library(BiocParallel)
register(MulticoreParam(4)) # Use 8 cores

# read Etv2 peaks at D1
source('chipseq.r'); 
peak_file <- 'analysis/etv2_pioneer/results/MEF_Dox_D1_Etv2_summits.bed'
gr <- read.table(peak_file, header = FALSE, sep = '\t')
gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]), peak_id = gr[, 4], log10pvalue = -gr[, 5])
gr <- resize(gr, width = 2000, fix = 'center')  # extending centered at the summit
gr <- add.seqinfo(gr, genome = 'mm10')

# Add the feature indicating whether or not there is significant Etv2 motif
library(chromVARmotifs) # https://github.com/GreenleafLab/chromVARmotifs
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
motif.set <- 'homer_pwms'
data(list = motif.set, package = 'chromVARmotifs')
h <- 'Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer(0.967)'
se <- SummarizedExperiment(rowRanges = gr)
motif_ix <- matchMotifs(get(motif.set)[h], se, genome = 'mm10')
values(gr)$Etv2_motif <- assays(motif_ix)$motifMatches[, 1]
source('analysis/etv2_pioneer/lr.r'); gr <- add_epigenetic_features(gr)
gr_file <- 'analysis/etv2_pioneer/results/MEF_Dox_D1_Etv2_summits.rds'
saveRDS(gr, gr_file)

# ----------------------------------------------------------------------------
# [2019-05-02] Add a feature whether or not the Etv2 peak is accessible in MEF
# ----------------------------------------------------------------------------
peak_file <- 'analysis/etv2_pioneer/results/ATAC_MEF_NoDox_peaks.narrowPeak'
atac <- read.table(peak_file, header = FALSE, sep = '\t')
atac <- GRanges(seqnames = atac[, 1], range = IRanges(atac[, 2], atac[, 3]))
source('chipseq.r'); atac <- add.seqinfo(atac, genome = 'mm10')



X <- values(gr)$epigenetic_features
X <- log(X + 1)
library(irlba); V <- irlba(scale(X), nu = 10, nv = 1)$u
library(Rtsne); y <- Rtsne(V)$Y

feature_file <- 'analysis/etv2_pioneer/data/Etv2_MEF_features.rds'
saveRDS(gr, file = feature_file)


# ----------------------------------------------------------------------------
# [2019-04-17] Merge the ATAC-seq BAM files for ES/EB and MEF data
# ----------------------------------------------------------------------------
dataset <- 'dataset=Etv2ATAC_version=20190228a'
source('sra.r'); d <- read.dataset(dataset, touch = FALSE)
d <- transform(d, sra.run.dir = sra.run.dir(run), bam.file = sprintf('%s/%s.dedup.bam', sra.run.result.dir(run), run))
d <- transform(d, group2 = gsub('_rep\\d$', '', group))
sp <- split(d[, 'bam.file'], list(d[, 'group2']))
bam_files <- sprintf('analysis/etv2_pioneer/data/ATAC_%s.bam', names(sp))
d2 <- data.frame(group = names(sp), bam_file = bam_files)
rownames(d2) <- d2[, 'group']
d2 <- transform(d2, ga_file = gsub('.bam', '.rds', bam_file))	# GAlignments file for each BAM
d2 <- transform(d2, vplot_file = gsub('.bam', '_Etv2peaks_vplot.rds', bam_file))	# Vplot density file for each BAM
d2 <- transform(d2, vplot_motif_file = gsub('.bam', '_Etv2motif_vplot.rds', bam_file))	# Vplot density file for each BAM
d2 <- transform(d2, vplot_Etv2_D1_MEF_file = gsub('.bam', '_Etv2peaks_D1_MEF_vplot.rds', bam_file))	# Vplot density file for each BAM


library(BiocParallel)
register(MulticoreParam(4)) 
source('sequencing.r'); bplapply(1:length(sp), function(i) merge.bam.files(sp[[i]], bam_files[i]))


# ----------------------------------------------------------------------------
# [2019-05-02] Remove the mitochrondria reads
# Not necessary since the % of Mt reads is relatively low
# ----------------------------------------------------------------------------
# d2 <- transform(d2, bam_noMt_file = gsub('.bam', '_noMt.bam', bam_file))
# i <- 1
# source('chipseq.r'); remove.mitochondrial.reads(d2[i, 'bam_file'], d2[i, 'bam_noMt_file'])


# ----------------------------------------------------------------------------
# [2019-05-02] Remove the blacklist regions
# http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# [2019-05-03] Read the BAM files for the ATAC-seq data and store them as RDS files
# Only read the necessary fields for the downstream analysis
# It appears that only isize is essential
# [2019-05-07] Generate the replicate level GA file
# ----------------------------------------------------------------------------
library(GenomicAlignments)
bam_param <- ScanBamParam(
	flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE, isDuplicate = FALSE), 
	what = c('isize')
)

# sample level data
d2 <- transform(d2, ga_file = gsub('.bam', '.rds', bam_file))
for (i in 1:nrow(d2)){
	cat(sprintf('processing %d/%d\n', i, nrow(d2)))
	atac <- readGAlignmentPairs(d2[i, 'bam_file'], param = bam_param)
	saveRDS(atac, d2[i, 'ga_file'])
}

# the same for the replicate level data
source('analysis/etv2_pioneer/helper.r'); d <- read_ATAC_dataset(sample_level = FALSE)
for (i in 4:5){
  cat(sprintf('processing %d/%d\n', i, nrow(d)))
  atac <- readGAlignmentPairs(d[i, 'bam.file'], param = bam_param)
  saveRDS(atac, d[i, 'ga_file'])
}



# ----------------------------------------------------------------------------
# [2019-05-02] Look at the global density plot
# [2019-05-07] Generate the global-level density plot, which will be used to normalize
# each V-plot
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); d <- read_ATAC_dataset(sample_level = FALSE, touch = TRUE)
lapply(1:nrow(d), function(i){
	cat(sprintf('processing %d/%d\n', i, nrow(d)))
	atac <- readRDS(d[i, 'ga_file'])
	read1 <- GenomicAlignments::first(atac)
	x <- abs(elementMetadata(read1)$isize)
	x <- as.data.frame(table(factor(x, 1:2000)))
	colnames(x) <- c('fragment_size', 'counts')
	write.table(x, d[i, 'fragment_size_file'], sep = '\t', quote = FALSE, row.names = FALSE)
})

D <- do.call('rbind', lapply(isizes, function(x) table(factor(x, 1:2000))))
D <- as.matrix(Diagonal(x = 1 / rowSums(D)) %*% D)
rownames(D) <- d2[, 'group']
D_file <- 'analysis/etv2_pioneer/results/ATAC_global_density.rds'
saveRDS(D, D_file)

x <- data.frame(density = c(as.matrix(D)), insert_size = rep(1:2000, each = nrow(D)), condition = rep(d2[, 'group'], 2000))
x <- transform(x, condition = factor(condition, c("EB_NoDox_D25", "EB_Dox_D25", "EB_Dox_D25_Flk1pos", "MEF_NoDox", "MEF_Dox_D1", "MEF_Dox_D2", "MEF_Dox_D7", "MEF_Dox_D7_Flk1pos")))
library(ggplot2); ggplot(x, aes(x = insert_size, y = density)) + geom_line() + facet_grid(~ condition ) + scale_y_continuous(trans = "log10") + ylim(c(1e2, 1e5))



# ----------------------------------------------------------------------------
# [2019-05-02] Look at the ATAC-seq read density plot at each condition
# ----------------------------------------------------------------------------
library(ATACseqQC)
library(GenomicAlignments)
library(BiocParallel)
library(SummarizedExperiment)
source('analysis/etv2_pioneer/helper.r')
register(MulticoreParam(4)) # Use 8 cores
tf <- 'Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer(0.967)'
source('analysis/etv2_pioneer/helper.r'); gr <- get_motifs(tf = tf, exclude_promoter = TRUE, exclude_exons = TRUE)
gr <- resize(gr, width = 100, fix = 'center')  # extending centered at the summit

d2 <- transform(d2, ga_file = gsub('.bam', '.rds', bam_file))
isizes <- lapply(1:nrow(d2), function(i){
	cat(sprintf('processing %d/%d: %s\n', i, nrow(d2), d2[i, 'group']))
	atac <- readRDS(d2[i, 'ga_file'])
	atac <- GAlignmentPairs(shiftReads(first(atac)), shiftReads(second(atac)))
	atac <- subsetByOverlaps(atac, gr)
	read1 <- GenomicAlignments::first(atac)
	abs(elementMetadata(read1)$isize)
})

D <- do.call('rbind', lapply(isizes, function(x) table(factor(x, 1:2000))))
D <- as.matrix(Diagonal(x = 1 / rowSums(D)) %*% D)
rownames(D) <- d2[, 'group']
#D_file <- 'analysis/etv2_pioneer/results/ATAC_density_Etv2_motifs.rds'
D_file <- 'analysis/etv2_pioneer/results/ATAC_density_Etv2_motifs_adjusted.rds'
saveRDS(D, D_file)

D0_file <- 'analysis/etv2_pioneer/results/ATAC_global_density.rds'
D0 <- readRDS(D0_file)
DD <- log(D + 1e-10) - log(D0 + 1e-10)

ylim <- 1500
DD <- t(log2(t(D[c('MEF_Dox_D1', 'MEF_Dox_D2', 'MEF_Dox_D7', 'MEF_Dox_D7_Flk1pos'), ] + 1e-5) / (D['MEF_NoDox', ] + 1e-5)))
#DD <- t(log2(t(D[c('EB_Dox_D25', 'EB_Dox_D25_Flk1pos'), ] + 1e-5) / (D['EB_NoDox_D25', ] + 1e-5)))
image(t(as.matrix(DD[, 1:ylim])), col = colorpanel(100, low = 'blue', mid = 'black', high = 'red'), breaks = c(-10, seq(-0.5, 0.5, length.out = 99), 10), axes = FALSE)
abline(v = c(180, 247) / ylim, col = 'gold', lwd = 2, lty = 2)
abline(v = c(315, 437) / ylim, col = 'yellow', lwd = 2, lty = 2)
abline(v = 100 / ylim, col = 'green', lwd = 2, lty = 2)
axis(3, seq(0, ylim, by = 250) / ylim)


# ----------------------------------------------------------------------------
# [2019-05-02] The V-plot analysis of Etv2 ChIP-seq peaks
# ----------------------------------------------------------------------------
library(GenomicAlignments)
library(BiocParallel)
library(SummarizedExperiment)

source('analysis/etv2_pioneer/helper.r'); d <- read_ATAC_dataset(sample_level = FALSE)
par(mfrow = c(4, 4))
xlim <- c(-2000, 2000); ylim <- c(50, 500)
#lapply(which(d[, 'group2'] %in% c('MEF_NoDox', 'MEF_Dox_D1')), function(i){
#lapply(1:nrow(d), function(i){
lapply(c(2, 3, 8, 16), function(i){
	cat(sprintf('processing %s\n', d[i, 'group']))
	atac <- readRDS(d[i, 'ga_file'])
	rel <- vplot(atac, shift(gr, 0), upstream = abs(xlim[1]), downstream = abs(xlim[2]))
	smoothScatter(rel, nbin = 256, xlim = xlim, ylim = ylim, main = sprintf('%s\n%s', label, d[i, 'group']))
})


#	saveRDS(rel, file = d2[i, 'vplot_file'])
	saveRDS(rel, file = d2[i, 'vplot_Etv2_D1_MEF_file'])



Z <- lapply(is, function(i){
	cat(sprintf('processing %s\n', rownames(d2)[i]))
	rel <- readRDS(d2[i, 'vplot_Etv2_D1_MEF_file'])
	rel <- transform(rel, FragmentLength = factor(FragmentLength, ylim[1]:ylim[2]), distanceToBindingSite = factor(distanceToBindingSite, xlim[1]:xlim[2]))
	as(table(rel[, 1], rel[, 2]), 'matrix')
})
Z <- lapply(Z, function(z) image.smooth(z / sum(z))$z)
names(Z) <- rownames(d2)[is]


#ZZ <- log(Z[['EB_Dox_D25_Flk1pos']] + 1e-10) - log(Z[['EB_NoDox_D25']] + 1e-10)
ZZ <- log(Z[['MEF_Dox_D2']] + 1e-10) - log(Z[['MEF_NoDox']] + 1e-10)
#ZZ <- log(Z[['MEF_Dox_D7']] + 1e-10) - log(Z[['MEF_NoDox']] + 1e-10)
#ZZ <- log(Z[['MEF_Dox_D7_Flk1pos']] + 1e-10) - log(Z[['MEF_NoDox']] + 1e-10)
ZZ <- image.smooth(ZZ, theta = 2)$z
image(ZZ, col = colorpanel(100, low = 'blue', mid = 'black', high = 'red'), breaks = c(-10, seq(-0.4, 0.4, length.out = 99), 10))


image(ZZ2, col = colorpanel(100, low = 'blue', mid = 'black', high = 'red'), breaks = c(-1, seq(-0.2, 0.2, length.out = 99), 1))

#h <- 'MEF_NoDox'
h <- 'MEF_Dox_D1'
image(Z[[h]], col = colorpanel(100, low = 'white', high = 'blue'), breaks = c(0, seq(0, 1e-5, length.out = 99), 1), main = h)

lapply(d2[is, 'group'], function(h){
	rel <- readRDS(d2[h, 'vplot_Etv2_D1_MEF_file'])
	smoothScatter(rel, nbin = 512, xlim = c(-500, 500), ylim = c(50, 500), main = h)
})

smoothScatter(rel, xlim = c(-500, 500), ylim = c(50, 500), colramp = Lab.palette)

#rel2 <- rel[rel$FragmentLength <= 500 & abs(rel$distanceToStart) < 500 , ]

h <- 'MEF_Dox_D7'
atac <- readRDS(d2[h, 'ga_file'])
rel <- vplot(atac, shift(gr, 100), upstream = 500, downstream = 500)
smoothScatter(rel, nbin = 256, xlim = c(-500, 500), ylim = c(50, 500), main = h)



Z <- as(table(rel2$FragmentLength, rel2$distanceToBindingSite), 'matrix')
image(t(Z), col = colorpanel(100, low = 'black', mid = 'green', high = 'yellow'))


library(MASS); kd <- kde2d(rel2[, 1], rel2[, 2], n = 100, lims = c(-200, 200, 70, 500))
image(kd, col = colorpanel(100, low = 'black', mid = 'green', high = 'yellow'))


atac <- subsetByOverlaps(atac, gr)	# ATAC-seq peaks overlapping with the motifs
source('analysis/etv2_pioneer/helper.r'); atac <- GAlignmentPairs(shiftReads(first(atac)), shiftReads(second(atac)))

x <- as(atac, 'GRanges')
w <- width(x)
x <- resize(x, width = 1, fix = 'center')
ol <- findOverlaps(x, gr)
bam.query <- bamIn[queryHits(ol)]
mt.ext.subject <- mt.ext[subjectHits(ol)]
rel <- getRelationship(bam.query, mt.ext.subject)



# ----------------------------------------------------------------------------
# [2019-05-06] Read ATAC-seq counts of all candidate Etv2 binding sites
# This will be used to generate the V-plot that focued on the Etv2 binding sites
# ----------------------------------------------------------------------------
library(GenomicAlignments)
library(BiocParallel)
library(SummarizedExperiment)
library(chromVAR)
register(MulticoreParam(2)) # Use 8 cores
tf <- 'Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer(0.967)'
source('analysis/etv2_pioneer/helper.r'); gr <- get_motifs(tf = tf, exclude_promoter = TRUE, exclude_exons = TRUE)
library(ChIPpeakAnno); gr <- reCenterPeaks(gr, width = 100)

# sampl level
source('analysis/etv2_pioneer/helper.r'); d2 <- read_ATAC_dataset(sample_level = TRUE)
source('chipseq.r'); se <- get_counts(d2[, 'bam_file'], gr)
se_file <- 'analysis/etv2_pioneer/data/se_ATAC_counts_on_Etv2_motifs.rds'
colData(se) <- DataFrame(d2)
saveRDS(se, se_file)

# repeat level
source('analysis/etv2_pioneer/helper.r'); d <- read_ATAC_dataset(sample_level = FALSE)
set.seed(1); source('chipseq.r'); se <- get_counts(d[, 'bam.file'], gr)
se_file <- 'analysis/etv2_pioneer/data/se_ATAC_counts_on_Etv2_motifs_replicate.rds'
colData(se) <- DataFrame(d)
saveRDS(se, se_file)


# ----------------------------------------------------------------------------
# [2019-05-06] Find the Etv2 motifs that meet the following criteria:
# * High at D1, low in MEF
# ----------------------------------------------------------------------------
library(GenomicAlignments)
library(BiocParallel)
library(SummarizedExperiment)
library(chromVAR)
register(MulticoreParam(2)) # Use 8 cores
se_file <- 'analysis/etv2_pioneer/data/se_ATAC_counts_on_Etv2_motifs_replicate.rds'	# ATAC-seq counts in each 100-bp Etv2 motifs
se <- readRDS(se_file)

h3k27ac <- read.table('analysis/etv2_pioneer/results/MEF_NoDox_H3K27ac_peaks.broadPeak', header = FALSE, sep = '\t')
h3k27ac <- GRanges(seqnames = h3k27ac[, 1], range = IRanges(h3k27ac[, 2], h3k27ac[, 3]))
non_enhancer <- granges(se) %outside% h3k27ac 

#library(DESeq2); dds <- DESeqDataSet(se, design = ~ group2)
#dds <- estimateSizeFactors(dds)
#dds <- DESeq(dds)
dds_file <- 'analysis/etv2_pioneer/data/se_ATAC_counts_on_Etv2_motifs_replicate_dds.rds'	# the DESeq2 object
#saveRDS(dds, dds_file)
dds <- readRDS(dds_file)
library(DESeq2); res <-  as.data.frame(results(dds, contrast = c('group2', 'MEF_Dox_D1', 'MEF_NoDox')))
X <- assays(se)$counts; colnames(X) <- colData(se)$group2
#n <- res[, 'log2FoldChange'] > 2 & res[, 'baseMean'] < 5
#n <- res[, 'log2FoldChange'] > 1 & res[, 'baseMean'] < 5
#n <- res[, 'log2FoldChange'] > 1 & res[, 'baseMean'] < 2
#n <- res[, 'baseMean'] < 2
#n <- res[, 'log2FoldChange'] > 1 & res[, 'baseMean'] < 2 & non_enhancer; label <- 'fc>1 & base<2 & non_enhancer'
#n <- res[, 'log2FoldChange'] > 2 & res[, 'baseMean'] < 2 & non_enhancer; label <- 'fc>2 & base<2 & non_enhancer'
n <- res[, 'log2FoldChange'] > 3 & res[, 'baseMean'] < 2 & non_enhancer; label <- 'fc>3 & base<2 & non_enhancer'
gr <- granges(se)[which(n), ]



library(ChIPpeakAnno); etv2 <- reCenterPeaks(etv2, width = 2000)
gr <- subsetByOverlaps(gr, etv2)	# ATAC-seq peaks overlapping with the motifs
source('analysis/etv2_pioneer/helper.r'); d <- read_ATAC_dataset(sample_level = FALSE)
source('chipseq.r'); se <- get_counts(d[, 'bam.file'], gr, mc.cores = 4)
colData(se) <- DataFrame(d)
se_file <- 'analysis/etv2_pioneer/data/se_ATAC_counts_on_Etv2_motifs.rds'
saveRDS(se, se_file)


treatment <- 'MEF_Dox_D2'
control <- 'MEF_NoDox'
m <- colData(se)$group2 %in% c(treatment, control)
n <- rowSums(assays(se)$counts[, colData(se)$group2 %in% control] < 5) == 2
gr <- granges(se[n])

lapply(which(m), function(i){
	cat(sprintf('processing %s\n', d[i, 'group'])
	atac <- readRDS(d[i, 'ga_file'])


# ----------------------------------------------------------------------------
# [2019-05-10] Calling nucleosomes by Danpos
# The command should run at a "bioinfo" conda environment
# ----------------------------------------------------------------------------
dataset <- 'dataset=Chronis_version=20190424a'
source('sra.r'); d <- read.dataset(dataset, touch = TRUE)
source('sra.r'); d <- transform(d, bam.file = sra.run.alignment.bam.file(run))

# convert the MNase-seq's BigWig to wig
# DANPOS accept wig as their input
#command <- sprintf('bigWigToWig /panfs/roc/scratch/gongx030/ncbi/sra/SRR50/SRR5077/SRR5077669/pileup.bw /panfs/roc/scratch/gongx030/ncbi/sra/SRR50/SRR5077/SRR5077669/pileupwig')
#command

#command <- sprintf('cp /panfs/roc/scratch/gongx030/ncbi/sra/SRR50/SRR5077/SRR5077669/pileup.wig %s/rlib/analysis/etv2_pioneer/data/Chronis/MNase', Sys.getenv('HOME'))

command <- sprintf('bedtools bamtobed -i /panfs/roc/scratch/gongx030/ncbi/sra/SRR50/SRR5077/SRR5077669/SRR5077669.bam > /panfs/roc/scratch/gongx030/ncbi/sra/SRR50/SRR5077/SRR5077669/SRR5077669.bed')
command

command <- sprintf('cp /panfs/roc/scratch/gongx030/ncbi/sra/SRR50/SRR5077/SRR5077669/SRR5077669.bed %s/rlib/analysis/etv2_pioneer/data/Chronis/MNase', Sys.getenv('HOME'))
command

input_dir <- sprintf('%s/rlib/analysis/etv2_pioneer/data/Chronis/MNase', Sys.getenv('HOME'))
bed_file <- sprintf('%s/SRR5077669.bed', input_dir)
dir.create(input_dir)
#command <- sprintf('bedtools bamtobed -i %s > %s', bam.file, bed_file)

#system(sprintf('cp %s %s', bam.file, input_dir))
command <- sprintf('python $HOME/src/danpos-2.2.2/danpos.py dpos %s', input_dir)
command


# ----------------------------------------------------------------------------
# [2019-05-12] Read the results from Danpos
# [2019-05-12] read_excel from readxl package cannot recognize the xls file
# [2019-05-12] The position.xls file are not in a XLS format; it is a plain tsv file
# ----------------------------------------------------------------------------
position_file <- sprintf('%s/rlib/analysis/etv2_pioneer/data/Chronis/MNase/result/pooled/home_garrydj_gongx030_rlib_analysis_etv2_pioneer_data_Chronis_MNase.smooth.positions.xls', Sys.getenv('HOME'))
#library(readxl)
#x <- read_excel(position_file)
x <- read.table(position_file, header = TRUE, sep = '\t')
source('chipseq.r'); gr <- GRanges(seqnames = x[, 1], range = IRanges(x[, 2], x[, 3]), smt_pos = x[, 4], smt_value = x[, 5], fuzziness_score = x[, 6])
gr <- add.seqinfo(gr, genome = 'mm10')
saveRDS(gr, file = 'analysis/etv2_pioneer/data/Chronis/MNase.smooth.positions.rds')


# ----------------------------------------------------------------------------
# [2019-05-12] Convert the .wig output from Danpos to the bigwig format
# ----------------------------------------------------------------------------
wig_file <- sprintf('%s/rlib/analysis/etv2_pioneer/data/Chronis/MNase/result/pooled/home_garrydj_gongx030_rlib_analysis_etv2_pioneer_data_Chronis_MNase.smooth.wig', Sys.getenv('HOME'))
bw_file <- sprintf('%s/rlib/analysis/etv2_pioneer/data/Chronis/MNase.smooth.bw', Sys.getenv('HOME'))
source('chipseq.r'); command <- sprintf('wigToBigWig %s %s %s', wig_file, genome.file('mm10', extend = TRUE), bw_file)
system(command)


# ----------------------------------------------------------------------------
# [2019-05-12] Use a nucleosome centric approach to visualize the Etv2 ChIP-seq peaks
# [2019-05-12] store the normalizeToMatrix for each condition to save time
# [2019-05-23] Use normalizeToMatrix_batch() from helper.r to generate the intermeidate N2M files
# [2019-05-12] Draw the heatmap for Etv2 peaks at D1 as well as other histone codes
# [2019-05-12] Use EnrichedHeatmap for drawing heatmap
# [2019-05-23] Cluster the data based on the histone code
# [2019-05-27]  Remove the central dotted line
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); etv2 <- read_Etv2_peaks_in_D1_MEF(exclude_promoter = TRUE, exclude_exons = FALSE)
nucleosome <- readRDS('analysis/etv2_pioneer/data/Chronis/MNase.smooth.positions.rds')
nucleosome <- GRanges(seqnames = seqnames(nucleosome), range = IRanges(mcols(nucleosome)$smt_pos, mcols(nucleosome)$smt_pos))
# the nucleosome summit with 200bp of Etv2 submmit
peaks <- subsetByOverlaps(nucleosome, reCenterPeaks(etv2, width = 200))
peaks <- add.seqinfo(peaks, genome = 'mm10')

# cluster the intervals based on histone modification peaks
gs <- c('MEF_H3K27ac', 'MEF_Brg1')
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
bw_files <- bw_files[gs]
bed_files <- rep(NA, length(bw_files)); names(bed_files) <- gs
i <- gs %in% c('MEF_Brg1', 'MEF_P300')
bed_files[i] <- gsub('_treat_pileup.bw', '_peaks.narrowPeak', bw_files[i])
i <- gs %in% c('MEF_H3K9me3', 'MEF_H3K27me3', 'MEF_H3K27ac')
bed_files[i] <- gsub('_treat_pileup.bw', '_peaks.broadPeak', bw_files[i])

S <- do.call('cbind', lapply(1:length(gs), function(i){
	gr <- read.table(bed_files[i], header = FALSE, sep = '\t')
	gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]))
	gr <- add.seqinfo(gr, genome = 'mm10')
	1:length(peaks) %in% as.matrix(findOverlaps(reCenterPeaks(peaks, width = 100 * 2), gr))[, 1]
}))
class(S) <- 'numeric'; colnames(S) <- gs
code <- Reduce('paste0', lapply(1:ncol(S), function(i) S[, i]))
#set.seed(1); km <- kmeans(S, 4)$cluster
set.seed(1); km <- S[, 'MEF_H3K27ac'] > 0


extend <- 1000; w <- 50
gs2 <- c('MEF_Dox_D1_Etv2', 'MEF_nucleosome', 'MEF_Brg1', 'MEF_H3K9me3', 'MEF_H3K27me3', 'MEF_MNase', 'MEF_H3K27ac', 'MEF_P300', 'MEF_H3', 'MEF_ATAC')
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
source('analysis/etv2_pioneer/helper.r'); n2m_files <- normalizeToMatrix_batch(peaks, peak_set = 'MEF_nucleosome_near_Etv2_peaks', bw_files[gs2], extend = extend, w = w, mc.cores = 4)


mat <- lapply(n2m_files, readRDS); names(mat) <- names(n2m_files)
col_fun <- lapply(mat, function(m) colorRamp2(quantile(m, c(0, 0.99)), c('white', 'blue')))
names(col_fun) <- names(n2m_files)

#i <- sample.int(nrow(mat[[1]]), 5000)
i <- 1:nrow(mat[[1]])
ta <- HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = 1:2, lty = 1), yaxis_facing = 'left'))
axis_name <- c('-1k', 'summit', '+1k')
h <- EnrichedHeatmap(mat[['MEF_H3']][i, ], split = factor(km[i], c(FALSE, TRUE)), col = col_fun[['MEF_H3']], name = 'H3', top_annotation = ta, axis_name = axis_name, pos_line = FALSE) + 
EnrichedHeatmap(mat[['MEF_nucleosome']][i, ], col = col_fun[['MEF_nucleosome']], name = 'DANPOS', top_annotation = ta, axis_name = axis_name, pos_line = FALSE) + 
EnrichedHeatmap(mat[['MEF_MNase']][i, ], col = col_fun[['MEF_MNase']], name = 'Mnase', top_annotation = ta, axis_name = axis_name, pos_line = FALSE) + 
EnrichedHeatmap(mat[['MEF_H3K27ac']][i, ], col = col_fun[['MEF_H3K27ac']], name = 'H3K27ac', top_annotation = ta, axis_name = axis_name, pos_line = FALSE) + 
EnrichedHeatmap(mat[['MEF_H3K9me3']][i, ], col = col_fun[['MEF_H3K9me3']], name = 'H3K9me3', top_annotation = ta, axis_name = axis_name, pos_line = FALSE) + 
EnrichedHeatmap(mat[['MEF_ATAC']][i, ], col = col_fun[['MEF_ATAC']], name = 'ATAC', top_annotation = ta, axis_name = axis_name, pos_line = FALSE) + 
EnrichedHeatmap(mat[['MEF_Dox_D1_Etv2']][i, ], col = colorRamp2(quantile(mat[['MEF_Dox_D1_Etv2']], c(0, 0.99)), c('white', 'red')), name = 'Etv2', top_annotation = ta, axis_name = axis_name, pos_line = FALSE)
draw(h, heatmap_legend_side = 'right')


# ----------------------------------------------------------------------------
# [2019-05-29] Address Ken's question
# I hate to ask you to do more work, but just so we can compare things directly with what we have done, could you please try centering on ETV sites but rank ordered by MNase tags +/- 100 base pairs from the centers?  
# That is what we have found to reliably predict in vitro nucleosome binding.
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); peaks <- read_Etv2_peaks_in_D1_MEF(exclude_promoter = TRUE, exclude_exons = FALSE)
peaks <- reCenterPeaks(peaks, width = 1)

gs <- c('MEF_H3K27ac')
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
bw_files <- bw_files[gs]
bed_files <- gsub('_treat_pileup.bw', '_peaks.broadPeak', bw_files)
gr <- read.table(bed_files, header = FALSE, sep = '\t')
gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]))
gr <- add.seqinfo(gr, genome = 'mm10')
with_H3K27ac <- 1:length(peaks) %in% as.matrix(findOverlaps(reCenterPeaks(peaks, width = 500 * 2), gr))[, 1]


extend <- 1000; w <- 50
gs <- c('MEF_H3', 'MEF_Dox_D1_Etv2', 'MEF_nucleosome', 'MEF_MNase', 'MEF_MNase2', 'MEF_H3K4me2', 'MEF_H3K4me1', 'MEF_H3K27ac', 'MEF_ATAC')
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
source('analysis/etv2_pioneer/helper.r'); n2m_files <- normalizeToMatrix_batch(peaks, peak_set = 'Etv2_peaks_D1', bw_files[gs], extend = extend, w = w, mc.cores = 4)

mat <- lapply(n2m_files, readRDS); names(mat) <- names(n2m_files)
col_fun <- lapply(mat, function(m) colorRamp2(quantile(m, c(0, 0.99)), c('white', 'blue')))
names(col_fun) <- names(n2m_files)

#i <- sample.int(nrow(mat[[1]]), 5000)
#i <- (!with_H3K27ac)
#i <- which(!with_H3K27ac)
#order_by <- order(rowSums(mat[['MEF_nucleosome']][!with_H3K27ac, (20 - 2):(20 + 2)]), decreasing = TRUE)
#i <- sample.int(nrow(mat[[1]]), 5000)
i <- 1:nrow(mat[[1]])
order_by <- order(rowSums(mat[['MEF_nucleosome']][i, (20 - 2):(20 + 2)]), decreasing = TRUE)

axis_name <- c('-1k', 'summit', '+1k')
h <- EnrichedHeatmap(mat[['MEF_Dox_D1_Etv2']][i, ], row_order = order_by, col = colorRamp2(quantile(mat[['MEF_Dox_D1_Etv2']], c(0, 0.99)), c('white', 'red')), name = 'Etv2', axis_name = axis_name, pos_line = FALSE) +
EnrichedHeatmap(mat[['MEF_H3']][i, ], col = col_fun[['MEF_H3']], name = 'H3', axis_name = axis_name, pos_line = FALSE) +
EnrichedHeatmap(mat[['MEF_nucleosome']][i, ], col = col_fun[['MEF_nucleosome']], name = 'DANPOS', axis_name = axis_name, pos_line = FALSE) +
EnrichedHeatmap(mat[['MEF_MNase']][i, ], col = col_fun[['MEF_MNase']], name = 'Mnase', axis_name = axis_name, pos_line = FALSE) + 
EnrichedHeatmap(mat[['MEF_MNase2']][i, ], col = col_fun[['MEF_MNase2']], name = 'Mnase2', axis_name = axis_name, pos_line = FALSE) + 
EnrichedHeatmap(mat[['MEF_H3K27ac']][i, ], col = col_fun[['MEF_H3K27ac']], name = 'H3K27ac', axis_name = axis_name, pos_line = FALSE)
#EnrichedHeatmap(mat[['MEF_ATAC']][i, ], col = col_fun[['MEF_ATAC']], name = 'ATAC', axis_name = axis_name, pos_line = FALSE) + 
#EnrichedHeatmap(mat[['MEF_H3K4me2']][i, ], col = col_fun[['MEF_H3K4me2']], name = 'H3K4me2', axis_name = axis_name, pos_line = FALSE) + 
#EnrichedHeatmap(mat[['MEF_H3K4me1']][i, ], col = col_fun[['MEF_H3K4me1']], name = 'H3K4me1', axis_name = axis_name, pos_line = FALSE) 
draw(h, heatmap_legend_side = 'right')


# ----------------------------------------------------------------------------
# [2019-05-29] Compare the MNase-seq and H3 ChIP-seq in MEF
# Address Ken's concern: It's interesting that the MNase data looks like it's more enriched in the lower half of the "low" group, opposite of H3.
# ----------------------------------------------------------------------------
nucleosome <- readRDS('analysis/etv2_pioneer/data/Chronis/MNase.smooth.positions.rds')
nucleosome <- GRanges(seqnames = seqnames(nucleosome), range = IRanges(mcols(nucleosome)$smt_pos, mcols(nucleosome)$smt_pos))
# the nucleosome summit with 200bp of Etv2 submmit
nucleosome <- add.seqinfo(nucleosome, genome = 'mm10')
set.seed(1); peaks <- nucleosome[sample.int(length(nucleosome), 50000)]

extend <- 1000; w <- 50
gs2 <- c('MEF_nucleosome', 'MEF_MNase', 'MEF_H3')
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
source('analysis/etv2_pioneer/helper.r'); n2m_files <- normalizeToMatrix_batch(peaks, peak_set = 'random_nucleosome_MEF', bw_files[gs2], extend = extend, w = w, mc.cores = 4)

mat <- lapply(n2m_files, readRDS); names(mat) <- names(n2m_files)
col_fun <- lapply(mat, function(m) colorRamp2(quantile(m, c(0, 0.99)), c('white', 'blue')))
names(col_fun) <- names(n2m_files)

i <- sample.int(nrow(mat[[1]]), 5000)
#i <- 1:nrow(mat[[1]])

axis_name <- c('-1k', 'summit', '+1k')
h <- EnrichedHeatmap(mat[['MEF_MNase']][i, ], col = col_fun[['MEF_MNase']], name = 'MNase', axis_name = axis_name, pos_line = FALSE) +
EnrichedHeatmap(mat[['MEF_nucleosome']][i, ], col = col_fun[['MEF_nucleosome']], name = 'DANPOS', axis_name = axis_name, pos_line = FALSE)
EnrichedHeatmap(mat[['MEF_H3']][i, ], col = col_fun[['MEF_H3']], name = 'H3', axis_name = axis_name, pos_line = FALSE) +
draw(h, heatmap_legend_side = 'right')






# ----------------------------------------------------------------------------
# [2019-05-14] Add Choi's Etv2 peaks to S3
# [2019-05-14] The bigbed (bb) format is not supported by rtracklayer.  Use bed format instead
# ----------------------------------------------------------------------------
#v5_file <- 'analysis/single_cell/etv2_chipseq_V5_mm10.bb'
v5_file <- 'analysis/single_cell/etv2_chipseq_V5_mm10.sorted.bed'
source('s3.r');  s3.backup(v5_file, 's3://etv2_pioneer', make_public = TRUE)

#pab_file <- 'analysis/single_cell/etv2_chipseq_polyab_mm10.bb'
pab_file <- 'analysis/single_cell/etv2_chipseq_polyab_mm10.sorted.bed'
source('s3.r');  s3.backup(pab_file, 's3://etv2_pioneer', make_public = TRUE)


# ----------------------------------------------------------------------------
# [2019-05-06] Compare the Etv2 MEF ChIP-seq Choi's Etv2 ChIP-seq data
# [2019-05-14] Use read_Etv2_peaks to read Etv2 peaks from MEF and ES/EB; Use findOverlapsOfPeaks
# from ChIPpeakAnno to merge peaks
# ----------------------------------------------------------------------------
library(ChIPpeakAnno)
library(EnrichedHeatmap)
library(circlize)
source('analysis/etv2_pioneer/helper.r'); peaks <- read_Etv2_peaks(width = 500)
peaks <- reCenterPeaks(peaks, width = 1)

#extend <- 1000; w <- 10; mc.cores <- 2
extend <- 1000; w <- 50; mc.cores <- 4
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
n2m_files <- sprintf('analysis/etv2_pioneer/results/normalizeToMatrix/peaks=Etv2_group=%s_extend=%d_w=%d.rds', names(bw_files), extend, w)
names(n2m_files) <- names(bw_files)

# store the normalizeToMatrix for each condition to save time
gs <- names(bw_files)
library(parallel); mclapply(gs, function(g){
	if (!file.exists(n2m_files[g])){
		cat(sprintf('processing %s\n', g))
		cvg <- import(bw_files[g], which = trim(reduce(reCenterPeaks(peaks, width = extend * 2))))	# returned as a GRanges object
		mat <- normalizeToMatrix(cvg, peaks, extend = extend, value_column = 'score', mean_mode = 'w0', w = w)
		cat(sprintf('writing %s\n', n2m_files[g]))
		saveRDS(mat, n2m_files[g])
	}
}, mc.cores = mc.cores)



# ----------------------------------------------------------------------------
# [2019-05-15] Draw the heatmap for Etv2 peaks at MEF D1/D2/D7 as well as Choi's ES/EB
# [2019-05-27] Focus on the H3K27ac-low Etv2 peaks at D1 MEF
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); peaks <- read_Etv2_peaks_in_D1_MEF(exclude_promoter = FALSE, exclude_exons = FALSE, MEF_H3K27ac = FALSE)
peaks <- reCenterPeaks(peaks, width = 1)
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()

gs <- c('MEF_Dox_D1_Etv2', 'MEF_H3K27ac', 'MEF_Dox_D2_Etv2', 'MEF_Dox_D7_Etv2', 'MEF_D1_H3K27ac', 'MEF_D2_H3K27ac', 'MEF_D7_H3K27ac', 'EB_D4_Etv2_pAb', 'EB_D4_Etv2_V5', 'MEF_ATAC', 'MEF_NoDox_ATAC', 'MEF_D1_ATAC', 'MEF_D2_ATAC', 'MEF_D7_ATAC', 'MEF_D7_Flk1pos_ATAC', 'EB_NoDox_D25', 'EB_Dox_D25', 'EB_Dox_D25_Flk1pos')
extend <- 1000; w <- 50
source('analysis/etv2_pioneer/helper.r'); n2m_files <- normalizeToMatrix_batch(peaks, peak_set = 'Etv2_peaks_D1', bw_files[gs], extend = extend, w = w, mc.cores = 4)


mat <- lapply(n2m_files, readRDS); names(mat) <- names(n2m_files)
col_fun <- lapply(mat, function(m) colorRamp2(quantile(m, c(0, 0.99)), c('white', 'blue')))
names(col_fun) <- names(n2m_files)

set.seed(1); i <- sample.int(nrow(mat[[1]]), 5000)
#i <- 1:nrow(mat[[1]])
#ta <- HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = 1:4, lty = 1), yaxis_facing = 'left'))
#axis_name <- c('-1k', 'summit', '+1k')
h <- EnrichedHeatmap(mat[['MEF_Dox_D1_Etv2']][i, ], col = col_fun[['MEF_Dox_D1_Etv2']], name = 'Etv2_D1') + 
EnrichedHeatmap(mat[['MEF_Dox_D2_Etv2']][i, ], col = col_fun[['MEF_Dox_D2_Etv2']], name = 'Etv2_D2') + 
EnrichedHeatmap(mat[['MEF_Dox_D7_Etv2']][i, ], col = col_fun[['MEF_Dox_D7_Etv2']], name = 'Etv2_D7') + 
EnrichedHeatmap(mat[['MEF_ATAC']][i, ], col = col_fun[['MEF_ATAC']], name = 'MEF_ATAC') + 
EnrichedHeatmap(mat[['MEF_NoDox_ATAC']][i, ], col = col_fun[['MEF_NoDox_ATAC']], name = 'ATAC_D0') + 
EnrichedHeatmap(mat[['MEF_D1_ATAC']][i, ], col = col_fun[['MEF_D1_ATAC']], name = 'ATAC_D1') + 
EnrichedHeatmap(mat[['MEF_D2_ATAC']][i, ], col = col_fun[['MEF_D2_ATAC']], name = 'ATAC_D2') + 
EnrichedHeatmap(mat[['MEF_D7_ATAC']][i, ], col = col_fun[['MEF_D7_ATAC']], name = 'ATAC_D7') 

#EnrichedHeatmap(mat[['MEF_H3K27ac']][i, ], col = col_fun[['MEF_H3K27ac']], name = 'MEF_H3K27ac') + 
#EnrichedHeatmap(mat[['EB_D4_Etv2_pAb']][i, ], col = col_fun[['EB_D4_Etv2_pAb']], name = 'Etv2_pAb') + 
#EnrichedHeatmap(mat[['EB_D4_Etv2_V5']][i, ], col = col_fun[['EB_D4_Etv2_V5']], name = 'Etv2_V5')

draw(h, heatmap_legend_side = 'right')




#extend <- 1000; w <- 50; gs <- c('MEF_Dox_D1_Etv2', 'MEF_Dox_D2_Etv2', 'MEF_Dox_D7_Etv2', 'EB_D4_Etv2_V5')
#extend <- 1000; w <- 50; gs <- c('MEF_Dox_D1_Etv2', 'MEF_Dox_D2_Etv2', 'MEF_Dox_D7_Etv2', 'EB_D4_Etv2_V5', 'MEF_ATAC', 'MEF_Brg1', 'MEF_H3', 'MEF_H3K27ac', 'MEF_P300', 'MEF_H3K9me3', 'MEF_H3K27me3', 'MEF_H3K36me3', 'MEF_H3K9ac', 'MEF_H3K79me2', 'MEF_H3K4me2', 'MEF_H3K4me1')
#extend <- 1000; w <- 50; gs <- c('MEF_Dox_D1_Etv2', 'MEF_Dox_D2_Etv2', 'MEF_Dox_D7_Etv2', 'EB_D4_Etv2_V5', 'MEF_H3K27ac', 'MEF_NoDox_H3K27ac', 'MEF_D1_H3K27ac', 'MEF_D2_H3K27ac', 'MEF_D7_H3K27ac')
#extend <- 1000; w <- 50; gs <- c('MEF_ATAC', 'MEF_NoDox_ATAC', 'MEF_D1_ATAC', 'MEF_D2_ATAC', 'MEF_D7_ATAC', 'MEF_D7_Flk1pos_ATAC', 'EB_NoDox_D25', 'EB_Dox_D25', 'EB_Dox_D25_Flk1pos')
extend <- 1000; w <- 50; gs <- c('MEF_ATAC', 'MEF_NoDox_rep1', 'MEF_NoDox_rep2', 'MEF_Dox_D1_rep1', 'MEF_Dox_D1_rep2', 'MEF_Dox_D2_rep1', 'MEF_Dox_D2_rep2', 'MEF_Dox_D7_rep1', 'MEF_Dox_D7_rep2', 'MEF_Dox_D7_Flk1pos_rep1', 'MEF_Dox_D7_Flk1pos_rep2', 'EB_NoDox_D25_rep1', 'EB_NoDox_D25_rep2', 'EB_Dox_D25_rep1', 'EB_Dox_D25_rep2', 'EB_Dox_D25_Flk1pos_rep1', 'EB_Dox_D25_Flk1pos_rep2')
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
n2m_files <- sprintf('analysis/etv2_pioneer/results/normalizeToMatrix/peaks=Etv2_group=%s_extend=%d_w=%d.rds', names(bw_files), extend, w)
names(n2m_files) <- names(bw_files)
mat <- lapply(gs, function(g) readRDS(n2m_files[g])); names(mat) <- gs
source('analysis/etv2_pioneer/helper.r'); cols <- get_bigwig_color()


set.seed(1); i <- sample.int(nrow(X), 10000)
#i <- 1:nrow(X)
#i <- mcols(peaks)$cluster %in% c('0001', '1001', '0101', '0011', '1101', '1011', '0111', '1111')
#i <- mcols(peaks)$cluster %in% c('1000', '0100', '0010', '0001', '1001', '1101', '1011', '1111')
#i <- mcols(peaks)$cluster %in% c('1000', '0100')
j <- rowSums(do.call('cbind', lapply(names(mat), function(g) rowSums(mat[[g]]) > quantile(rowSums(mat[[g]]), 0.999)))) == 0
col_fun <- lapply(names(mat), function(g) colorRamp2(quantile(mat[[g]][i & j, ], c(0, 0.99)), c('black', cols[g]))); names(col_fun) <- names(mat)
h1 <- Heatmap(X[i & j, ], name = 'group',  col = c('white', 'black'), cluster_rows = FALSE, cluster_columns = FALSE, split = cls[i & j], width = unit(30, 'mm'))
h_list <- lapply(gs, function(g) EnrichedHeatmap(mat[[g]][i & j, ], col = col_fun[[g]], name = g))
draw(h1 + Reduce('+', h_list), heatmap_legend_side = 'bottom')



# ----------------------------------------------------------------------------
# [2019-05-12] Compar the distribution of peaks at genomic locations
# ----------------------------------------------------------------------------
library(ChIPpeakAnno)
library(EnrichedHeatmap)
library(circlize)
source('analysis/etv2_pioneer/helper.r'); peaks <- read_Etv2_peaks(width = 500)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

d <- lapply(1:length(metadata(peaks)$ol$peaklist), function(i){
	cat(sprintf('%d/%d\n', i, length(metadata(peaks)$ol$peaklist)))
	assignChromosomeRegion(metadata(peaks)$ol$peaklist[[i]], nucleotideLevel = FALSE, precedence = c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Exons", "Introns"), TxDb= TxDb.Mmusculus.UCSC.mm10.knownGene)
})

X <- as(metadata(peaks)$ol$venn_cnt, 'matrix')[-1, 1:4]
y <- as(metadata(peaks)$ol$venn_cnt, 'matrix')[-1, 5]
i <- y > 500
Z <- do.call('rbind', lapply(d, function(xx) xx$percentage))

par(mar = c(5, 5, 5, 5))
peak_counts <- HeatmapAnnotation(
	number = anno_barplot(
		y[i], 
		which = 'row',
		bar_width = 1, 
		gp = gpar(col = "white", fill = "black"), 
		border = TRUE,
		axis_side = 'bottom',
		axis = TRUE,
		axis_param = list(at = c(0, 20000, 40000), labels = c('0', '20k', '40k')), 
	),
	which = 'row',
	width = unit(3, 'cm'),
	show_annotation_name = FALSE,
	name = 'time'
)
Heatmap(X[i, ], show_row_dend = FALSE, cluster_columns = FALSE, col = c('white', 'black')) + peak_counts + Heatmap(Z[i, ], show_column_dend = FALSE, name = 'genome')


# ----------------------------------------------------------------------------
# [2019-05-15] Generating repeat level ATAC-seq pileup (MACS2)
# There appears something wrong with the sample level data
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); d <- read_ATAC_dataset(sample_level = FALSE)
treatment_files <- d[, 'bam.file']
base.name <- sprintf('analysis/etv2_pioneer/results/ATAC_%s', d[, 'group'])

source('chipseq.r');
library(parallel); mclapply(1:nrow(d), function(i){
	 macs2.callpeak(treatment_files[i], base.name[i], format = 'BAMPE', genome = 'mm10', broad = FALSE, qvalue.cutoff = 0.05, fold.change = FALSE, update = TRUE, call.summits = TRUE, shift = -100, extsize = 200)
}, mc.cores = 2)

pileup_file <- sprintf('%s_treat_pileup.bw', base.name)
source('s3.r'); s3.backup(pileup_file, 's3://etv2_pioneer/',  make_public = TRUE)



# ----------------------------------------------------------------------------
# [2019-05-16] Generating repeat level H3K27ac pileup (MACS2)
# [2019-05-17] Put the MACS2 results to the $PROJECT_DIR/macs2 to save storage space
# $PROJECT_DIR is defined in helper.r
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); d <- read_H3K27ac_dataset(sample_level = FALSE)

pairs <- list(
	 list(name = 'EB_H3K27ac_12h_rep1', treatment = c("EB_H3K27ac_12h_rep1"), control = c('EB_Dox_12h_input')),
	 list(name = 'EB_H3K27ac_12h_rep2', treatment = c("EB_H3K27ac_12h_rep2"), control = c('EB_Dox_12h_input')),
	 list(name = 'EB_H3K27ac_3h_rep1', treatment = c("EB_H3K27ac_3h_rep1"), control = c('EB_Dox_3h_input')),
	 list(name = 'EB_H3K27ac_3h_rep2', treatment = c("EB_H3K27ac_3h_rep2"), control = c('EB_Dox_3h_input')),
	 list(name = 'EB_NoDox_12h_H3K27ac_rep1', treatment = c("EB_NoDox_12h_H3K27ac_rep1"), control = c('EB_Dox_12h_input')),
	 list(name = 'EB_NoDox_12h_H3K27ac_rep2', treatment = c("EB_NoDox_12h_H3K27ac_rep2"), control = c('EB_Dox_12h_input')),
	 list(name = 'EB_NoDox_3h_H3K27ac_rep1', treatment = c("EB_NoDox_3h_H3K27ac_rep1"), control = c('EB_Dox_3h_input')),
	 list(name = 'EB_NoDox_3h_H3K27ac_rep2', treatment = c("EB_NoDox_3h_H3K27ac_rep2"), control = c('EB_Dox_3h_input')),
	 list(name = 'MEF_Dox_D1_H3K27ac_rep1', treatment = c("MEF_Dox_D1_H3K27ac_rep1"), control = c('MEF_Dox_D1_input_rep1')),
	 list(name = 'MEF_Dox_D1_H3K27ac_rep2', treatment = c("MEF_Dox_D1_H3K27ac_rep2"), control = c('MEF_Dox_D1_input_rep1')),
	 list(name = 'MEF_Dox_D2_H3K27ac_rep1', treatment = c("MEF_Dox_D2_H3K27ac_rep1"), control = c('MEF_Dox_D2_input_rep1')),
	 list(name = 'MEF_Dox_D2_H3K27ac_rep2', treatment = c("MEF_Dox_D2_H3K27ac_rep2"), control = c('MEF_Dox_D2_input_rep1')),
	 list(name = 'MEF_Dox_D7_H3K27ac_rep1', treatment = c("MEF_Dox_D7_H3K27ac_rep1"), control = c('MEF_Dox_D7_input_rep1')),
	 list(name = 'MEF_Dox_D7_H3K27ac_rep2', treatment = c("MEF_Dox_D7_H3K27ac_rep2"), control = c('MEF_Dox_D7_input_rep1')),
	 list(name = 'MEF_NoDox_H3K27ac_rep1', treatment = c("MEF_NoDox_H3K27ac_rep1"), control = c('MEF_NoDox_input_rep1')),
	 list(name = 'MEF_NoDox_H3K27ac_rep2', treatment = c("MEF_NoDox_H3K27ac_rep2"), control = c('MEF_NoDox_input_rep1'))
)

source('chipseq.r');
library(parallel); mclapply(1:length(pairs), function(i){
	base.name <- sprintf('%s/macs2/%s', PROJECT_DIR, pairs[[i]]$name)
	treatment_files <- d[d$group %in% pairs[[i]]$treatment, 'bam.file']
	control_files <- d[d$group %in% pairs[[i]]$control, 'bam.file']
	macs2.callpeak(treatment_files, base.name, control_files, format = 'BAM', genome = 'mm10', broad = TRUE, qvalue.cutoff = 0.05, fold.change = FALSE, update = TRUE, call.summits = FALSE, shift = 0, extsize = 200, keep.dup = 'all')
}, mc.cores = 4)

pileup_file <- sprintf('%s_treat_pileup.bw', base.name)
source('s3.r'); s3.backup(pileup_file, 's3://etv2_pioneer/',  make_public = TRUE)


# ----------------------------------------------------------------------------
# [2019-06-02] Generating sampl3 level H3K27ac pileup (MACS2)
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); d <- read_H3K27ac_dataset(sample_level = TRUE)

treatment_files <- d[d$group2 %in% c('EB_H3K27ac_12h'), 'bam.file']
control_files <- d[d$group2 %in% c('EB_Dox_12h_input'), 'bam.file']
base.name <- sprintf('%s/macs2/%s', PROJECT_DIR, 'EB_H3K27ac_12h')
macs2.callpeak(treatment_files, base.name, control_files, format = 'BAMPE', genome = 'mm10', broad = TRUE, qvalue.cutoff = 0.05, fold.change = FALSE, update = TRUE, call.summits = FALSE, shift = 0, extsize = 200, keep.dup = 'all')



# ----------------------------------------------------------------------------
# [2019-05-16] Read ATAC-seq peaks
# [2019-05-16] Generate the N2M intermediate results
# [2019-05-16] Try different cutoffs for ATAC-seq peaks.  The default q<0.05 cutoff generate many relative faint peaks
# Try to use a more stringent cutoff (e.g. q < 1e-10) to get a cleaner heatmap
# [2019-05-21] Generate the N2M for H3K27ac data
# ----------------------------------------------------------------------------
library(ChIPpeakAnno)
library(EnrichedHeatmap)
library(circlize)
source('analysis/etv2_pioneer/helper.r'); peaks <- read_ATAC_peaks(width = 500, log10qvalue = 0)
peaks <- reCenterPeaks(peaks, width = 1)

extend <- 1000; w <- 50; mc.cores <- 2
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
n2m_files <- sprintf('analysis/etv2_pioneer/results/normalizeToMatrix/peaks=ATAC_group=%s_extend=%d_w=%d.rds', names(bw_files), extend, w)	# use q < 0.05 as the cutoff (default)
#n2m_files <- sprintf('analysis/etv2_pioneer/results/normalizeToMatrix/peaks=ATACq10_group=%s_extend=%d_w=%d.rds', names(bw_files), extend, w)	# use q<1e-10 as the cutoff
#n2m_files <- sprintf('analysis/etv2_pioneer/results/normalizeToMatrix/peaks=ATACq5_group=%s_extend=%d_w=%d.rds', names(bw_files), extend, w)	# use q<1e-10 as the cutoff
names(n2m_files) <- names(bw_files)

# store the normalizeToMatrix for each condition to save time
#gs <- names(bw_files)
#gs <- c('MEF_NoDox_ATAC', 'MEF_D1_ATAC', 'MEF_D2_ATAC', 'MEF_D7_ATAC', 'MEF_D7_Flk1pos_ATAC', 'EB_NoDox_D25', 'EB_Dox_D25', 'EB_Dox_D25_Flk1pos')
#gs <- c('MEF_NoDox_H3K27ac', 'MEF_D2_H3K27ac', 'MEF_D7_H3K27ac')
#gs <- c('MEF_Dox_D2_Etv2', 'MEF_Dox_D7_Etv2')
#gs <- c('EB_D4_Etv2_V5')
gs <- c('MEF_Dox_D1_Etv2')
library(parallel); mclapply(gs, function(g){
	if (!file.exists(n2m_files[g])){
		cat(sprintf('processing %s\n', g))
		cvg <- import(bw_files[g], which = trim(reduce(reCenterPeaks(peaks, width = extend * 2))))  # returned as a GRanges object
		mat <- normalizeToMatrix(cvg, peaks, extend = extend, value_column = 'score', mean_mode = 'w0', w = w)
		cat(sprintf('writing %s\n', n2m_files[g]))
		saveRDS(mat, n2m_files[g])
	}
}, mc.cores = mc.cores)


# ----------------------------------------------------------------------------
# [2019-05-16] Generate the N2M intermediate results
# Accessibility change during the MEF induction
# ----------------------------------------------------------------------------
library(EnrichedHeatmap)
library(circlize)
source('analysis/etv2_pioneer/helper.r'); peaks <- read_ATAC_peaks(width = 500, log10qvalue = 0)
peaks <- reCenterPeaks(peaks, width = 1)
X <- mcols(peaks)$group; class(X) <- 'numeric'


#extend <- 1000; w <- 50; gs <- c('MEF_NoDox_ATAC', 'MEF_D2_ATAC', 'MEF_D7_Flk1pos_ATAC')
#extend <- 1000; w <- 50; gs <- c('MEF_NoDox_ATAC', 'MEF_D1_ATAC', 'MEF_D2_ATAC', 'MEF_D7_ATAC', 'MEF_D7_Flk1pos_ATAC', 'EB_NoDox_D25', 'EB_Dox_D25', 'EB_Dox_D25_Flk1pos', 'MEF_Dox_D1_Etv2', 'MEF_Dox_D2_Etv2', 'MEF_Dox_D7_Etv2', 'EB_D4_Etv2_V5')
#extend <- 1000; w <- 50; gs <- c('MEF_NoDox_ATAC', 'MEF_NoDox_H3K27ac', 'MEF_D2_ATAC', 'MEF_D2_H3K27ac', 'MEF_D7_ATAC', 'MEF_D7_H3K27ac')
#extend <- 1000; w <- 50; gs <- c('MEF_NoDox_ATAC', 'MEF_D2_ATAC', 'MEF_D7_Flk1pos_ATAC', 'EB_Dox_D25_Flk1pos')
#extend <- 1000; w <- 50; gs <- c('MEF_NoDox_ATAC', 'MEF_D1_ATAC', 'EB_NoDox_D25', 'EB_Dox_D25')
#extend <- 1000; w <- 50; gs <- c('MEF_Dox_D1_Etv2', 'MEF_Dox_D2_Etv2', 'MEF_Dox_D7_Etv2', 'MEF_NoDox_ATAC', 'MEF_D1_ATAC', 'MEF_D2_ATAC', 'MEF_D7_ATAC', 'MEF_D7_Flk1pos_ATAC')
#extend <- 1000; w <- 50; gs <- c('MEF_NoDox_ATAC', 'MEF_D1_ATAC', 'MEF_Dox_D1_Etv2', 'MEF_D2_ATAC', 'MEF_Dox_D2_Etv2', 'MEF_D7_ATAC', 'MEF_Dox_D7_Etv2', 'MEF_D7_Flk1pos_ATAC', 'EB_Dox_D25_Flk1pos', 'EB_D4_Etv2_V5')
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
n2m_files <- sprintf('analysis/etv2_pioneer/results/normalizeToMatrix/peaks=ATAC_group=%s_extend=%d_w=%d.rds', names(bw_files), extend, w)
#n2m_files <- sprintf('analysis/etv2_pioneer/results/normalizeToMatrix/peaks=ATACq10_group=%s_extend=%d_w=%d.rds', names(bw_files), extend, w)
names(n2m_files) <- names(bw_files)
mat <- lapply(gs, function(g) readRDS(n2m_files[g])); names(mat) <- gs
source('analysis/etv2_pioneer/helper.r'); cols <- get_bigwig_color()


#set.seed(1); cls <- kmeans(X, 6)$cluster
#set.seed(1); i <- 1:nrow(X) %in% sample.int(nrow(X), 10000)
#ss <- c('MEF_NoDox_ATAC', 'MEF_D2_ATAC', 'MEF_D7_Flk1pos_ATAC')
#ss <- c('MEF_NoDox_ATAC', 'MEF_D2_ATAC', 'MEF_D7_Flk1pos_ATAC')
ss <- c('MEF_NoDox_ATAC', 'MEF_D2_ATAC', 'MEF_D7_Flk1pos_ATAC')
#set.seed(1); i <- rowSums(X[, ss]) > 0 & rowSums(X[, ss]) < length(ss) 
j <- rowSums(do.call('cbind', lapply(names(mat), function(g) rowSums(mat[[g]]) > quantile(rowSums(mat[[g]]), 0.999)))) == 0
set.seed(1); cls <- kmeans(X[i & j, ss], 6)$cluster
col_fun <- lapply(names(mat), function(g) colorRamp2(quantile(mat[[g]][i & j, ], c(0, 0.999)), c('black', cols[g]))); names(col_fun) <- names(mat)
h2 <- Heatmap(X[i & j, ss], name = 'group',  col = c('white', 'black'), width = unit(50, 'mm'), cluster_columns = FALSE)
h1 <- EnrichedHeatmap(mat[[1]][i & j, ], col = col_fun[[1]], name = gs[1], split = cls)
h_list <- lapply(2:length(gs), function(h) EnrichedHeatmap(mat[[h]][i & j, ], col = col_fun[[h]], name = gs[h]))
draw(h1 + Reduce('+', h_list) + h2, heatmap_legend_side = 'bottom')




i <- X[, 'MEF_NoDox_ATAC'] == 0 & X[, 'MEF_D1_ATAC'] == 0 & X[, 'MEF_D2_ATAC'] == 0 & (X[, 'MEF_D7_ATAC'] == 1  | X[, 'MEF_D7_Flk1pos_ATAC'] == 1)
j <- rowSums(do.call('cbind', lapply(names(mat), function(g) rowSums(mat[[g]]) > quantile(rowSums(mat[[g]]), 0.999)))) == 0
col_fun <- lapply(names(mat), function(g) colorRamp2(quantile(mat[[g]][i & j, ], c(0, 0.999)), c('black', cols[g]))); names(col_fun) <- names(mat)
h1 <- EnrichedHeatmap(mat[[1]][i & j, ], col = col_fun[[1]], name = gs[1])
h_list <- lapply(2:length(gs), function(h) EnrichedHeatmap(mat[[h]][i & j, ], col = col_fun[[h]], name = gs[h]))
draw(h1 + Reduce('+', h_list), heatmap_legend_side = 'bottom')






# show Etv2 help closing some clusters
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
n2m_files <- sprintf('analysis/etv2_pioneer/results/normalizeToMatrix/peaks=ATAC_group=%s_extend=%d_w=%d.rds', names(bw_files), extend, w)
names(n2m_files) <- names(bw_files)
mat <- lapply(gs, function(g) readRDS(n2m_files[g])); names(mat) <- gs
source('analysis/etv2_pioneer/helper.r'); cols <- get_bigwig_color()


cls_number <- 6
ss <- c('MEF_NoDox_ATAC', 'MEF_D2_ATAC', 'MEF_D7_Flk1pos_ATAC')
set.seed(1); i <- rowSums(X[, ss]) > 0 & rowSums(X[, ss]) < length(ss)
j <- rowSums(do.call('cbind', lapply(names(mat), function(g) rowSums(mat[[g]]) > quantile(rowSums(mat[[g]]), 0.999)))) == 0
set.seed(1); cls <- kmeans(X[i & j, ss], 6)$cluster

j <- rowSums(do.call('cbind', lapply(names(mat), function(g) rowSums(mat[[g]]) > quantile(rowSums(mat[[g]]), 0.999)))) == 0
set.seed(1); cls <- kmeans(X[i & j, ss], 6)$cluster
col_fun <- lapply(names(mat), function(g) colorRamp2(quantile(mat[[g]][i & j, ], c(0, 0.999)), c('black', cols[g]))); names(col_fun) <- names(mat)

h1 <- EnrichedHeatmap(mat[[1]][i & j, ], col = col_fun[[1]], name = gs[1], split = cls)

h_list <- lapply(2:length(gs), function(h) EnrichedHeatmap(mat[[h]][i & j, ], col = col_fun[[h]], name = gs[h]))
draw(h1 + Reduce('+', h_list), heatmap_legend_side = 'bottom')


# ----------------------------------------------------------------------------
# [2019-05-17] chromVAR validation of the sample level ATAC-seq data
# [2019-05-21] The results are similar to the repeat-level data
# It should be OK to go ahead using the sample level data for the downstream analysis
# ----------------------------------------------------------------------------
library(chromVAR)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
source('analysis/etv2_pioneer/helper.r'); peaks <- read_ATAC_peaks(width = 500, log10qvalue = 10)
source('analysis/etv2_pioneer/helper.r'); d <- read_ATAC_dataset(sample_level = TRUE)
source('chipseq.r'); se <- get_counts(d[, 'bam_file'], peaks)	# this step takes ~ 1 hour
se_file <- 'analysis/etv2_pioneer/results/se_ATAC_repeat_level.rds'
saveRDS(se, se_file)

se_file <- 'analysis/etv2_pioneer/results/se_ATAC_repeat_level.rds'
se <- readRDS(se_file)

peaks <- reCenterPeaks(peaks, width = 500)
motif.set <- 'homer_pwms'
data(list = motif.set, package = 'chromVARmotifs')
se <- addGCBias(se, genome = BSgenome.Mmusculus.UCSC.mm10)
motif_ix <- matchMotifs(get(motif.set), se, genome = 'mm10')
dev <- computeDeviations(object = se, annotations = motif_ix)
v <- computeVariability(dev)

h <- v[, 'p_value_adj'] < 0.05
Z <- assays(dev)$z[h, ]
pca <- prcomp(Z)

library(wordcloud); textplot(pca$rotation[, 1], pca$rotation[, 2], d[, 'group'], xpd = TRUE)



# ----------------------------------------------------------------------------
# [2019-05-17] Look at the motifs near the Etv2 occupied nucleosome
# [2019-05-22] Only consider the motifs overlapping with the ATAC footprinting
# [2019-05-26] Consider two groups of Etv2-nearby-nucleosomes: the H3K27ac low and high region
# The TF distance relationship might be different for these two groups of Etv2 peaks
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r')
#width <- 600; MEF_H3K27ac <- NULL; motif_file <- sprintf('%s/X_motif.rds', PROJECT_DIR)
width <- 600; MEF_H3K27ac <- TRUE; motif_file <- sprintf('%s/X_motif_H3K27ac=TRUE.rds', PROJECT_DIR)
#width <- 600; MEF_H3K27ac <- FALSE; motif_file <- sprintf('%s/X_motif_H3K27ac=FALSE', PROJECT_DIR)

source('analysis/etv2_pioneer/helper.r'); etv2 <- read_Etv2_peaks_in_D1_MEF(exclude_promoter = FALSE, exclude_exons = FALSE, MEF_H3K27ac = MEF_H3K27ac)
nucleosome <- readRDS('analysis/etv2_pioneer/data/Chronis/MNase.smooth.positions.rds')
nucleosome <- GRanges(seqnames = seqnames(nucleosome), range = IRanges(mcols(nucleosome)$smt_pos, mcols(nucleosome)$smt_pos))
peaks <- subsetByOverlaps(nucleosome, reCenterPeaks(etv2, width = 200))
peaks <- add.seqinfo(peaks, genome = 'mm10')
source('analysis/etv2_pioneer/helper.r'); pwms <- get_PWM_in_MEF()

# The motifs near the nucleosome surrounding the Etv2 peaks at D1 MEF
motif_ix <- matchMotifs(pwms, reCenterPeaks(peaks, width = width), genome = 'mm10', out = 'positions')
hint <- read.table('analysis/etv2_pioneer/data/ATAC_MEF_Dox_D1_HINT.bed', header = FALSE, sep = '\t')
hint <- GRanges(seqnames = hint[, 1], ranges = IRanges(hint[, 2], hint[, 3]))

library(parallel); X_motif <- do.call('rbind', mclapply(1:length(motif_ix), function(i){
	cat(sprintf('%d/%d\n', i, length(motif_ix)))
	gr <- subsetByOverlaps(motif_ix[[i]], hint)
	gr <- add.seqinfo(gr, genome = 'mm10')
	gr <- reCenterPeaks(gr, 50)
	colMeans(as(coverage(gr)[reCenterPeaks(peaks, width)], 'matrix'))
}, mc.cores = 4))
saveRDS(X_motif, motif_file)



# Look at the nucleosomes NOT near any Etv2 peaks
#peaks_neg <- nucleosome[!nucleosome %in% reCenterPeaks(etv2, width = 2000)]	# nucleosome outside of Etv2 peaks at D1 MEF
source('analysis/etv2_pioneer/helper.r'); etv2_motif <- get_motifs(tf = 'Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer(0.967)', exclude_promoter = FALSE, exclude_exons = FALSE)
peaks_neg <- subsetByOverlaps(nucleosome, reCenterPeaks(etv2_motif, width = 200))	# nucleosome near the Etv2 motifs
set.seed(1); peaks_neg <- peaks_neg[sample.int(length(peaks_neg), 100000)]	# randomly sample 100k
peaks_neg <- add.seqinfo(peaks_neg, genome = 'mm10')
motif_ix_neg <- matchMotifs(pwms, reCenterPeaks(peaks_neg, width = width), genome = 'mm10', out = 'positions')

library(parallel); X_motif_neg <- do.call('rbind', mclapply(1:length(motif_ix_neg), function(i){
	cat(sprintf('%d/%d\n', i, length(motif_ix_neg)))
	gr <- subsetByOverlaps(motif_ix_neg[[i]], hint)
	gr <- add.seqinfo(gr, genome = 'mm10')
	gr <- reCenterPeaks(gr, 50)
	colMeans(as(coverage(gr)[reCenterPeaks(peaks_neg, width)], 'matrix'))
}, mc.cores = 4))
saveRDS(X_motif_neg, sprintf('%s/X_motif_neg.rds', PROJECT_DIR))



# ----------------------------------------------------------------------------
# [2019-05-22] Visualize the positions of motifs to the nucleosome center
# These motifs must overlap with the ATAC-seq footprinting at D1 MEF
# ----------------------------------------------------------------------------
#width <- 600; MEF_H3K27ac <- TRUE; motif_file <- sprintf('%s/X_motif_H3K27ac=TRUE.rds', PROJECT_DIR)
width <- 600; MEF_H3K27ac <- FALSE; motif_file <- sprintf('%s/X_motif_H3K27ac=FALSE', PROJECT_DIR)
X_motif <- readRDS(motif_file)
rownames(X_motif) <- names(pwms)
center <- 300; extend <- 250
source('analysis/etv2_pioneer/helper.r'); Yp <- matrix2normalizedMatrix(X_motif, center = center, extend = extend)

y <- apply(Yp, 1, which.max) - 74
i <- order(y)
gs <- c('Sp1', 'Pitx2', 'Id1', 'Tgif1', 'Smad4', 'Jun', 'Nfia', 'Arid5b', 'Tfe3', 'Gata6', 'Gata2', 'Rela', 'E2f1', 'Klf4', 'Myc', 'Fosl1', 'Jun', 'Ctcf', 'Etv2', 'Stat3', 'Hoxb2', 'Runx1', 'Twist1', 'Meis1', 'Ebf1', 'Yy1', 'Sox11', 'Foxp1', 'Mef2c', 'Mbd1', 'Prrx2', 'Irf2', 'Tbp', 'Foxk1', 'Six4', 'Hmg20b')
at <- which(rownames(Yp) %in% gs)
h <- EnrichedHeatmap(Yp, row_order = rev(i), show_row_names = FALSE, col = colorRamp2(quantile(Yp, c(0, 0.5, 0.99)), c('blue', 'white', 'red'))) +
rowAnnotation(link = row_anno_link(at = at, labels = rownames(Yp)[at]), width = unit(3, 'cm'))
draw(h, heatmap_legend_side = 'bottom')


# ----------------------------------------------------------------------------
# [2019-05-26] Jointly Visualize the positions of motifs to the nucleosome center 
# on H3K27ac present and absent regions
# ----------------------------------------------------------------------------
center <- 300; extend <- 250
width <- 600; MEF_H3K27ac <- TRUE; motif_file <- sprintf('%s/X_motif_H3K27ac=TRUE.rds', PROJECT_DIR)
X_motif <- readRDS(motif_file)
rownames(X_motif) <- names(pwms)
source('analysis/etv2_pioneer/helper.r'); Yp <- matrix2normalizedMatrix(X_motif, center = center, extend = extend)

width <- 600; MEF_H3K27ac <- FALSE; motif_file <- sprintf('%s/X_motif_H3K27ac=FALSE', PROJECT_DIR)
X_motif <- readRDS(motif_file)
rownames(X_motif) <- names(pwms)
source('analysis/etv2_pioneer/helper.r'); Yn <- matrix2normalizedMatrix(X_motif, center = center, extend = extend)

yn <- apply(Yn, 1, which.max)
gs <- c('Sp1', 'Pitx2', 'Id1', 'Tgif1', 'Smad4', 'Jun', 'Nfia', 'Arid5b', 'Tfe3', 'Gata6', 'Gata2', 'Rela', 'E2f1', 'Klf4', 'Myc', 'Fosl1', 'Jun', 'Ctcf', 'Etv2', 'Stat3', 'Hoxb2', 'Runx1', 'Twist1', 'Meis1', 'Ebf1', 'Yy1', 'Sox11', 'Foxp1', 'Mef2c', 'Mbd1', 'Prrx2', 'Irf2', 'Tbp', 'Foxk1', 'Six4', 'Hmg20b')
at <- which(rownames(Yn) %in% gs)
hn <- EnrichedHeatmap(Yn, row_order = rev(order(yn)), show_row_names = FALSE, col = colorRamp2(quantile(Yn, c(0, 0.5, 0.99)), c('blue', 'white', 'red')))
#hp <- EnrichedHeatmap(Yp, show_row_names = FALSE, col = colorRamp2(quantile(Yp, c(0, 0.5, 0.99)), c('blue', 'white', 'red')))
ra <- rowAnnotation(link = row_anno_link(at = at, labels = rownames(Yn)[at]), width = unit(3, 'cm'))
draw(hn + ra, heatmap_legend_side = 'bottom')


yp <- apply(Yp, 1, which.max)
d <- yn - yp
#i <- order(d)[1:10]
i <- order(d, decreasing = TRUE)[1:10]
x <- rbind(data.frame(tf = names(d)[i], H3K27ac = TRUE, distance = yp[i]), data.frame(tf = names(d)[i], H3K27ac = FALSE, distance = yn[i]))

library(ggplot2)
ggplot(x, aes(x = tf, y = distance, fill = H3K27ac)) + geom_bar(stat = 'identity', position=position_dodge()) + coord_flip() + theme(axis.text=element_text(size=16, face = 'bold')) + scale_fill_manual(values=c('black','red'))





# ----------------------------------------------------------------------------
# [2019-05-26] A heatmap of H3K27ac on/off nucleosomes
# with the density of Etv2 motif, Etv2 peak and Gata motif
# ----------------------------------------------------------------------------
g <- 'Gata2'
plot(as.numeric(Yp[g, ]), lwd = 2, col = 'red', ylim = range(c(Yp[g, ]), c(Yn[g, ])), type = 'l')
lines(as.numeric(Yn[g, ]), lwd = 2, col = 'black')
abline(v = 74, col = 'green', lwd = 2)







# ----------------------------------------------------------------------------
# [2019-05-21] Get genome-wise PWM hits across genome
# These results will be used for PIQ analysis
# e.g. Step 2 in https://bitbucket.org/thashim/piq-single/src/master/
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); pwms <- get_mouse_pwms_v2()
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Matrix)
# Looks like PIQ used all chromosomes from BSgenome.Mmusculus.UCSC.mm10, so no filtering is needed
# e.g. there are in total 66 chromosomes in BSgenome.Mmusculus.UCSC.mm10
#source('aux.r'); gr <- add.seqinfo(gr, genome = 'mm10')

pwm_file <- sprintf('%s/PIQ_pwms/mouse_pwms_v2.txt', PROJECT_DIR)
# WARNING: there is a blank line between every TFs, which will fail PIQ's pwmmatch.exact.r
# Need to manually remove all the blank lines for now
#library(BiocParallel)
#register(MulticoreParam(4)) # Use 8 cores
#source('analysis/etv2_pioneer/helper.r'); write_jaspar_for_PIQ(pwms, pwm_file)
mm_dir <- sprintf('%s/motif.matches', PROJECT_DIR)

lapply(9:length(pwms), function(i){
	tryCatch({
		cat(sprintf('%d/%d\n', i, length(pwms)))
		command <- sprintf('Rscript $PIQ_HOME/pwmmatch.exact.r $PIQ_HOME/common.r %s %d %s/', pwm_file, i, mm_dir)
		system(command)
	}, error = function(e) NULL)
})


# ----------------------------------------------------------------------------
# [2019-05-21]  Convert BAM to internal binary format (does not depend on choice of motif).
# https://bitbucket.org/thashim/piq-single/src/master/
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); d <- read_ATAC_dataset(sample_level = TRUE)
piq_bam_dir <- sprintf('%s/PIQ_bams', PROJECT_DIR)
d <- transform(d, piq_file = sprintf('%s/ATAC_%s.RData', piq_bam_dir, group))

library(parallel); mclapply(1:nrow(d), function(i){
	command <- sprintf('Rscript $PIQ_HOME/bam2rdata.r $PIQ_HOME/common.r %s %s', d[i, 'piq_file'], d[i, 'bam_file'])
	cat(sprintf('[%s] %s\n', Sys.time(), command))
	system(command)
}, mc.cores = 2)


# ----------------------------------------------------------------------------
# [2019-05-21] Calling TF binding sites using PIQ
# https://bitbucket.org/thashim/piq-single/src/master/
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); pwms <- get_mouse_pwms_v2()
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
source('analysis/etv2_pioneer/helper.r'); d <- read_ATAC_dataset(sample_level = TRUE)
piq_bam_dir <- sprintf('%s/PIQ_bams', PROJECT_DIR)
d <- transform(d, piq_file = sprintf('%s/ATAC_%s.RData', piq_bam_dir, group))
d <- transform(d, piq_pertf_dir = sprintf('%s/ATAC_%s_pertf', piq_bam_dir, group))


lapply(1:nrow(d), function(m){
	dir.create(d[m, 'piq_pertf_dir'])
	mclapply(1:length(pwms), function(tf){
		pwm_dir <- sprintf('%s/motif.matches', PROJECT_DIR)
		tmp <- tempdir()
		command <- sprintf('Rscript $PIQ_HOME/pertf.r $PIQ_HOME/common.r %s/ %s/ %s %s %s', pwm_dir, tmp, d[m, 'piq_pertf_dir'], d[m, 'piq_file'], tf)
		system(command)
	}, mc.cores = 2)
})


source('analysis/etv2_pioneer/helper.r'); pwms <- get_mouse_pwms_v2()
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BiocParallel)
register(MulticoreParam(4)) # Use 8 cores


# ----------------------------------------------------------------------------
# [2019-05-21] Read ATAC-seq peaks and save as a bed file for HINT
# ----------------------------------------------------------------------------
width <- 500; log10qvalue <- 0
source('analysis/etv2_pioneer/helper.r'); peaks <- read_ATAC_peaks(width = width, log10qvalue = log10qvalue)
peaks <- reCenterPeaks(peaks, width = width)
bed_file <- sprintf('analysis/etv2_pioneer/data/ATAC_peaks_log10qvalue=%d_width=%d.bed', log10qvalue, width)
write.table(as.data.frame(peaks)[, 1:3], bed_file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)


# ----------------------------------------------------------------------------
# [2019-05-21] Identifying the footprinting of in house ATAC-seq data using HINT
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); d <- read_ATAC_dataset(sample_level = TRUE)
width <- 500; log10qvalue <- 0
bed_file <- sprintf('analysis/etv2_pioneer/data/ATAC_peaks_log10qvalue=%d_width=%d.bed', log10qvalue, width)
d <- transform(d, hint_prefix = gsub('.bam', '_HINT', d[, 'bam_file']))

lapply(1:2, function(i){
	command <- sprintf('rgt-hint footprinting --paired-end --organism mm10 --atac-seq %s %s --output-prefix=%s', d[i, 'bam_file'], bed_file, d[i, 'hint_prefix'])
	cat(sprintf('[%s] %s\n', Sys.time(), command))
	system(command)
})

command <- sprintf('rgt-hint footprinting --paired-end --organism mm10 --atac-seq %s %s --output-prefix=%s', d[i, 'bam_file'], 'analysis/etv2_pioneer/data/test.bed', 'analysis/etv2_pioneer/data/test_HINT')



# ----------------------------------------------------------------------------
# [2019-05-17] Generate tracks for HINT's footprinting results
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); d <- read_ATAC_dataset(sample_level = TRUE)
d <- transform(d, hint_bed_file = gsub('.bam', '_HINT.bed', d[, 'bam_file']))
d <- transform(d, hint_bb_file = gsub('.bam', '_HINT.bb', d[, 'bam_file']))
source('filenames.r')
for (i in 1:nrow(d)){
	tfile <- tempfile()
	command <- sprintf('cat %s | cut -f 1,2,3,4 > %s', d[i, 'hint_bed_file'], tfile)
	cat(sprintf('[%s] %s\n', Sys.time(), command)); system(command)
	command <- sprintf('bedToBigBed  %s %s %s -as=$KENT_HOME/src/hg/lib/bed.as', tfile, genome.file('mm10', extend = TRUE), d[i, 'hint_bb_file'])
	cat(sprintf('[%s] %s\n', Sys.time(), command)); system(command)
}
source('s3.r');  s3.backup(d[, 'hint_bb_file'], 's3://etv2_pioneer/hint/', make_public = TRUE)


# ----------------------------------------------------------------------------
# [2019-05-17] Generate tracks for HINT's footprinting results
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
gs <- c('MEF_ATAC', 'MEF_H3', 'MEF_nucleosome', 'MEF_MNase', 'MEF_H3K27ac', 'MEF_H3K9me3')
bw_files <- bw_files[gs]
source('s3.r');  s3.backup(bw_files, 's3://etv2_pioneer/data/Chronis/', make_public = TRUE)



# ----------------------------------------------------------------------------
# [2019-05-27] A set of ATAC-seq peaks from MEF_D7_Flk1pos_ATAC and EB_Dox_D25_Flk1pos
# The peaks are grouped into three category 01, 10 and 11
# ----------------------------------------------------------------------------
width <- 200
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
gs <- c('MEF_NoDox_ATAC', 'MEF_D7_Flk1pos_ATAC', 'EB_Dox_D25_Flk1pos')
bw_files <- bw_files[gs]
bed_files <- gsub('_treat_pileup.bw', '_summits.bed', bw_files)
gr0 <- lapply(bed_files, function(file){
	gr <- read.table(file, header = FALSE, sep = '\t')
	gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]), log10qvalue = gr[, 5])
	gr <- reCenterPeaks(gr, width = width)
	gr <- add.seqinfo(gr, genome = 'mm10')
	gr <- unique(gr)
	gr
})
x <- findOverlapsOfPeaks(gr0[[1]], gr0[[2]])
peaks <- Reduce('c', x$peaklist)
peaks <- reCenterPeaks(peaks, width = 1)
group <- rep(1:length(x$peaklist), sapply(x$peaklist, length))

#gs <- c('MEF_Dox_D1_Etv2', 'MEF_H3K27ac', 'MEF_Dox_D2_Etv2', 'MEF_Dox_D7_Etv2', 'MEF_D1_H3K27ac', 'MEF_D2_H3K27ac', 'MEF_D7_H3K27ac', 'EB_D4_Etv2_pAb', 'EB_D4_Etv2_V5', 'MEF_ATAC', 'MEF_NoDox_ATAC', 'MEF_D1_ATAC', 'MEF_D2_ATAC', 'MEF_D7_ATAC', 'MEF_D7_Flk1pos_ATAC', 'EB_NoDox_D25', 'EB_Dox_D25', 'EB_Dox_D25_Flk1pos')
gs <- c('MEF_NoDox_ATAC', 'MEF_D1_ATAC', 'MEF_D2_ATAC', 'MEF_D7_ATAC', 'MEF_D7_Flk1pos_ATAC', 'EB_NoDox_D25', 'EB_Dox_D25', 'EB_Dox_D25_Flk1pos')
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
extend <- 1000; w <- 50
source('analysis/etv2_pioneer/helper.r'); n2m_files <- normalizeToMatrix_batch(peaks, peak_set = 'Flk1pos_ATAC', bw_files[gs], extend = extend, w = w, mc.cores = 4)

mat <- lapply(n2m_files, readRDS); names(mat) <- names(n2m_files)
col_fun <- lapply(mat, function(m) colorRamp2(quantile(m, c(0, 0.99)), c('white', 'blue')))
names(col_fun) <- names(n2m_files)

set.seed(1); i <- sample.int(nrow(mat[[1]]), 5000)
h <- EnrichedHeatmap(mat[['MEF_NoDox_ATAC']][i, ], split = factor(group[i], c(1, 2, 3)), col = col_fun[['MEF_NoDox_ATAC']], name = 'MEF_NoDox_ATAC') +
EnrichedHeatmap(mat[['MEF_D1_ATAC']][i, ], col = col_fun[['MEF_D1_ATAC']], name = 'MEF_D1_ATAC') +
EnrichedHeatmap(mat[['MEF_D2_ATAC']][i, ], col = col_fun[['MEF_D2_ATAC']], name = 'MEF_D2_ATAC') +
EnrichedHeatmap(mat[['MEF_D7_ATAC']][i, ], col = col_fun[['MEF_D7_ATAC']], name = 'MEF_D7_ATAC') +
EnrichedHeatmap(mat[['MEF_D7_Flk1pos_ATAC']][i, ], col = col_fun[['MEF_D7_Flk1pos_ATAC']], name = 'MEF_D7_Flk1pos_ATAC') +
EnrichedHeatmap(mat[['EB_NoDox_D25']][i, ], col = col_fun[['EB_NoDox_D25']], name = 'EB_NoDox_D25') +
EnrichedHeatmap(mat[['EB_Dox_D25']][i, ], col = col_fun[['EB_Dox_D25']], name = 'EB_Dox_D25') +
EnrichedHeatmap(mat[['EB_Dox_D25_Flk1pos']][i, ], col = col_fun[['EB_Dox_D25_Flk1pos']], name = 'EB_Dox_D25_Flk1pos')

draw(h, heatmap_legend_side = 'right')

# ----------------------------------------------------------------------------
# [2019-05-30] Testing the ideas that Etv2 targeting the nucleosome and the between region 
# have different motifs
# * Order the Etv2 D1 MEF peaks based on the MNase-seq signals
# * Looking at the enriched motifs at the botton and top quantile
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); peaks <- read_Etv2_peaks_in_D1_MEF(exclude_promoter = FALSE, exclude_exons = FALSE, MEF_H3K27ac = NULL)
peaks <- reCenterPeaks(peaks, width = 1)

extend <- 1000; w <- 50
gs <- c('MEF_H3', 'MEF_Dox_D1_Etv2', 'MEF_nucleosome', 'MEF_MNase')
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
source('analysis/etv2_pioneer/helper.r'); n2m_files <- normalizeToMatrix_batch(peaks, peak_set = 'All_Etv2_peaks_D1', bw_files[gs], extend = extend, w = w, mc.cores = 4)
mat <- lapply(n2m_files, readRDS); names(mat) <- names(n2m_files)

# Visualize the heatmap
col_fun <- lapply(mat, function(m) colorRamp2(quantile(m, c(0, 0.99)), c('white', 'blue')))
names(col_fun) <- names(n2m_files)
i <- 1:nrow(mat[[1]])
order_by <- order(rowSums(mat[['MEF_nucleosome']][i, (20 - 2):(20 + 2)]), decreasing = TRUE)
axis_name <- c('-1k', 'summit', '+1k')
h <- EnrichedHeatmap(mat[['MEF_Dox_D1_Etv2']][i, ], row_order = order_by, col = colorRamp2(quantile(mat[['MEF_Dox_D1_Etv2']], c(0, 0.99)), c('white', 'red')), name = 'Etv2', axis_name = axis_name, pos_line = FALSE) +
#EnrichedHeatmap(mat[['MEF_H3']][i, ], col = col_fun[['MEF_H3']], name = 'H3', axis_name = axis_name, pos_line = FALSE) +
EnrichedHeatmap(mat[['MEF_nucleosome']][i, ], col = col_fun[['MEF_nucleosome']], name = 'DANPOS', axis_name = axis_name, pos_line = FALSE) +
EnrichedHeatmap(mat[['MEF_MNase']][i, ], col = col_fun[['MEF_MNase']], name = 'Mnase', axis_name = axis_name, pos_line = FALSE)
draw(h, heatmap_legend_side = 'right')

# Find the motifs in Etv2 overlapped with MNase and not overlapped with the MNase peak
x <- rowSums(mat[['MEF_nucleosome']][, (20 - 2):(20 + 2)])
with_MNase <- x > quantile(x, 0.75)
with_MNase_peak_file <- sprintf('%s/macs2/Etv2_with_MNase_D1.bed', PROJECT_DIR)
write.table(as.data.frame(peaks)[with_MNase, 1:3], with_MNase_peak_file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

without_MNase <- x < quantile(x, 0.25)
without_MNase_peak_file <- sprintf('%s/macs2/Etv2_without_MNase_D1.bed', PROJECT_DIR)
write.table(as.data.frame(peaks)[without_MNase, 1:3], without_MNase_peak_file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)


# ----------------------------------------------------------------------------
# [2019-05-30] Use MEME to find the enriched motif in Etv2 peaks with or without MNase peaks
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r');
#peak_file <- sprintf('%s/macs2/Etv2_with_MNase_D1.bed', PROJECT_DIR)
peak_file <- 'analysis/single_cell/etv2_chipseq_V5_mm10.sorted.bed'
#peak_file <- sprintf('%s/macs2/Etv2_without_MNase_D1.bed', PROJECT_DIR)
fasta_file <- gsub('.bed', '.fa', peak_file)
output_dir <- gsub('.bed', '.meme', peak_file)
centrimo_output_dir <- gsub('.bed', '.centrimo', peak_file)
motif_file <- gsub('.bed', '.meme/dreme.txt', peak_file)

gr <- read.table(peak_file, header = FALSE, sep = '\t')
gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]))
gr <- reCenterPeaks(gr, width = 200)
s <- getSeq(BSgenome.Mmusculus.UCSC.mm10, gr)
names(s) <- 1:length(s)
writeXStringSet(s, fasta_file)
command <- sprintf('dreme -p %s -oc %s -png -mink 6 -maxk 7 -m 20', fasta_file, output_dir)
system(command)

command <- sprintf('centrimo %s %s --oc %s', fasta_file, motif_file, centrimo_output_dir)
system(command)

# Get the canonical Etv2 motif from Choi's data

# ----------------------------------------------------------------------------
# [2019-05-31] Make the sample level ATAC-seq data on S3
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); d <- read_ATAC_dataset(sample_level = TRUE)
source('s3.r'); s3.backup(d[, 'bam_file'], 's3://garrylab/bam/', make_public = TRUE)



# ----------------------------------------------------------------------------
# [2019-06-02] Look at the ATAC-seq signal change of the Etv2 binding in H3K27ac-low region
# in MEF
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); peaks <- read_Etv2_peaks_in_D1_MEF(exclude_promoter = TRUE, exclude_exons = FALSE, MEF_H3K27ac = FALSE)
peaks <- reCenterPeaks(peaks, width = 1)

extend <- 1000; w <- 50
gs <- c('MEF_NoDox_ATAC', 'MEF_D1_ATAC', 'MEF_D2_ATAC', 'MEF_D7_ATAC', 'MEF_D7_Flk1pos_ATAC', 'EB_Dox_D25_Flk1pos', 'MEF_D1_H3K27ac', 'MEF_D2_H3K27ac', 'MEF_D7_H3K27ac')
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
source('analysis/etv2_pioneer/helper.r'); n2m_files <- normalizeToMatrix_batch(peaks, peak_set = 'Etv2_peaks_D1_MEF_H3K27ac_low', bw_files[gs], extend = extend, w = w, mc.cores = 4)

mat <- lapply(n2m_files, readRDS); names(mat) <- names(n2m_files)
col_fun <- lapply(mat, function(m) colorRamp2(quantile(m, c(0, 0.99)), c('white', 'blue')))
names(col_fun) <- names(n2m_files)

i <- sample.int(nrow(mat[[1]]), 5000)
order_by <- order(rowSums(mat[['MEF_D7_H3K27ac']][i, (20 - 2):(20 + 2)]), decreasing = TRUE)

axis_name <- c('-1k', 'summit', '+1k')
h <- EnrichedHeatmap(mat[['MEF_NoDox_ATAC']][i, ], row_order = order_by, col = col_fun[['MEF_NoDox_ATAC']], name = 'MEF_NoDox_ATAC', axis_name = axis_name, pos_line = FALSE) +
EnrichedHeatmap(mat[['MEF_D1_ATAC']][i, ], col = col_fun[['MEF_D1_ATAC']], name = 'MEF_D1_ATAC', axis_name = axis_name, pos_line = FALSE) +
EnrichedHeatmap(mat[['MEF_D7_Flk1pos_ATAC']][i, ], col = col_fun[['MEF_D7_Flk1pos_ATAC']], name = 'MEF_D7_Flk1pos_ATAC', axis_name = axis_name, pos_line = FALSE) +
EnrichedHeatmap(mat[['EB_Dox_D25_Flk1pos']][i, ], col = col_fun[['EB_Dox_D25_Flk1pos']], name = 'EB_Dox_D25_Flk1pos', axis_name = axis_name, pos_line = FALSE) + 
EnrichedHeatmap(mat[['MEF_D1_H3K27ac']][i, ], col = col_fun[['MEF_D1_H3K27ac']], name = 'MEF_D1_H3K27ac', axis_name = axis_name, pos_line = FALSE) + 
EnrichedHeatmap(mat[['MEF_D2_H3K27ac']][i, ], col = col_fun[['MEF_D2_H3K27ac']], name = 'MEF_D2_H3K27ac', axis_name = axis_name, pos_line = FALSE) + 
EnrichedHeatmap(mat[['MEF_D7_H3K27ac']][i, ], col = col_fun[['MEF_D7_H3K27ac']], name = 'MEF_D7_H3K27ac', axis_name = axis_name, pos_line = FALSE) 
draw(h, heatmap_legend_side = 'right')


# ----------------------------------------------------------------------------
# [2019-06-02] Look at the ATAC-seq signal change between MEF_ATAC, MEF_D7_Flk1pos_ATAC and EB_Dox_D25_Flk1pos
# ----------------------------------------------------------------------------
source('analysis/etv2_pioneer/helper.r'); peaks <- read_ATAC_peaks(width = 200)
peaks <- reCenterPeaks(peaks, width = 1)

extend <- 1000; w <- 50
#gs <- c('MEF_NoDox_ATAC', 'MEF_D1_ATAC', 'MEF_D2_ATAC', 'MEF_D7_ATAC', 'MEF_D7_Flk1pos_ATAC', 'EB_NoDox_D25', 'EB_Dox_D25', 'EB_Dox_D25_Flk1pos', 'MEF_NoDox_H3K27ac')
#gs <- c('MEF_NoDox_ATAC', 'MEF_D7_ATAC', 'EB_Dox_D25', 'MEF_NoDox_H3K27ac', 'MEF_D7_H3K27ac', 'EB_H3K27ac_12h')
gs <- c('MEF_NoDox_ATAC', 'MEF_D7_Flk1pos_ATAC', 'EB_Dox_D25_Flk1pos', 'MEF_NoDox_H3K27ac', 'MEF_D7_H3K27ac', 'EB_H3K27ac_12h')
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
#source('analysis/etv2_pioneer/helper.r'); n2m_files <- normalizeToMatrix_batch(peaks, peak_set = 'ATAC_peaks_three_time_points', bw_files[gs], extend = extend, w = w, mc.cores = 4)
source('analysis/etv2_pioneer/helper.r'); n2m_files <- normalizeToMatrix_batch(peaks, peak_set = 'ATAC_peaks_three_time_points_Flk1pos', bw_files[gs], extend = extend, w = w, mc.cores = 4)

mat <- lapply(n2m_files, readRDS); names(mat) <- names(n2m_files)
cols <- c(
	'MEF_NoDox_ATAC' = 'blue', 
	'MEF_D7_ATAC' = 'blue', 
	'EB_Dox_D25' = 'blue', 
	'MEF_D7_Flk1pos_ATAC' = 'blue', 
	'EB_Dox_D25_Flk1pos' = 'blue', 
	'MEF_NoDox_H3K27ac' = 'red', 
	'MEF_D7_H3K27ac' = 'red', 
	'EB_H3K27ac_12h' = 'red'
)
col_fun <- lapply(names(mat), function(g) colorRamp2(quantile(mat[[g]], c(0, 0.98)), c('white', cols[g])))
names(col_fun) <- names(n2m_files)

i <- sample.int(nrow(mat[[1]]), 5000)
#i <- 1:nrow(mat[[1]])
split_by <- factor(mcols(peaks)$cluster[i], c('100', '010', '001', '110', '101', '011', '111'))

axis_name <- c('-1k', 'summit', '+1k')
h <- EnrichedHeatmap(mat[['MEF_NoDox_ATAC']][i, ], split = split_by, col = col_fun[['MEF_NoDox_ATAC']], name = 'MEF_NoDox_ATAC', axis_name = axis_name, pos_line = FALSE) +
#EnrichedHeatmap(mat[['MEF_D7_ATAC']][i, ], col = col_fun[['MEF_D7_ATAC']], name = 'MEF_D7_ATAC', axis_name = axis_name, pos_line = FALSE)  + 
#EnrichedHeatmap(mat[['EB_Dox_D25']][i, ], col = col_fun[['EB_Dox_D25']], name = 'EB_Dox_D25', axis_name = axis_name, pos_line = FALSE)  + 
EnrichedHeatmap(mat[['MEF_D7_Flk1pos_ATAC']][i, ], col = col_fun[['MEF_D7_Flk1pos_ATAC']], name = 'MEF_D7_Flk1pos_ATAC', axis_name = axis_name, pos_line = FALSE)  + 
EnrichedHeatmap(mat[['EB_Dox_D25_Flk1pos']][i, ], col = col_fun[['EB_Dox_D25_Flk1pos']], name = 'EB_Dox_D25_Flk1pos', axis_name = axis_name, pos_line = FALSE)  + 
EnrichedHeatmap(mat[['MEF_NoDox_H3K27ac']][i, ], col = col_fun[['MEF_NoDox_H3K27ac']], name = 'MEF_NoDox_H3K27ac', axis_name = axis_name, pos_line = FALSE) +
EnrichedHeatmap(mat[['MEF_D7_H3K27ac']][i, ], col = col_fun[['MEF_D7_H3K27ac']], name = 'MEF_D7_H3K27ac', axis_name = axis_name, pos_line = FALSE) +
EnrichedHeatmap(mat[['EB_H3K27ac_12h']][i, ], col = col_fun[['EB_H3K27ac_12h']], name = 'EB_H3K27ac_12h', axis_name = axis_name, pos_line = FALSE) 
draw(h, heatmap_legend_side = 'right')


i <- mcols(peaks)$cluster %in% '011'
order_by <- order(rowSums(mat[['MEF_NoDox_ATAC']][i, (20 - 2):(20 + 2)]), decreasing = TRUE)

axis_name <- c('-1k', 'summit', '+1k')
h <- EnrichedHeatmap(mat[['MEF_NoDox_ATAC']][i, ], row_order = order_by, col = col_fun[['MEF_NoDox_ATAC']], name = 'MEF_NoDox_ATAC', axis_name = axis_name, pos_line = FALSE) +
EnrichedHeatmap(mat[['MEF_D7_ATAC']][i, ], col = col_fun[['MEF_D7_ATAC']], name = 'MEF_D7_ATAC', axis_name = axis_name, pos_line = FALSE)  + 
EnrichedHeatmap(mat[['EB_Dox_D25']][i, ], col = col_fun[['EB_Dox_D25']], name = 'EB_Dox_D25', axis_name = axis_name, pos_line = FALSE)  + 
EnrichedHeatmap(mat[['MEF_NoDox_H3K27ac']][i, ], col = col_fun[['MEF_NoDox_H3K27ac']], name = 'MEF_NoDox_H3K27ac', axis_name = axis_name, pos_line = FALSE) +
EnrichedHeatmap(mat[['MEF_D7_H3K27ac']][i, ], col = col_fun[['MEF_D7_H3K27ac']], name = 'MEF_D7_H3K27ac', axis_name = axis_name, pos_line = FALSE) +
EnrichedHeatmap(mat[['EB_H3K27ac_12h']][i, ], col = col_fun[['EB_H3K27ac_12h']], name = 'EB_H3K27ac_12h', axis_name = axis_name, pos_line = FALSE) 
draw(h, heatmap_legend_side = 'right')


# ----------------------------------------------------------------------------
# [2019-08-22] Calling nucleosome on MEF ATAC-seq data
# ----------------------------------------------------------------------------




