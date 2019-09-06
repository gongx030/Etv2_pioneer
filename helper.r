library(ComplexHeatmap)
library(SummarizedExperiment)
library(ChIPpeakAnno)
library(futile.logger)
library(parallel)
library(EnrichedHeatmap)
library(circlize)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
source('chipseq.r')
source('sra.r'); 

PROJECT_DIR <- sprintf('%s/etv2_pioneer', Sys.getenv('TMPDIR'))

# touch everything in the project dir
# find /panfs/roc/scratch/gongx030/etv2_pioneer  -type f -exec touch {} +

flog.info(sprintf('PROJECT_DIR: %s', PROJECT_DIR))
flog.info(sprintf('MACS2 results: %s/macs2', PROJECT_DIR))

atacseq_group2bg <- function(group){
	g2bg <- c(
		"MEF_NoDox_rep1" = 'black', 
		"MEF_NoDox_rep2" = 'black',
		"MEF_Dox_D1_rep1" = 'blue',
		"MEF_Dox_D1_rep2" = 'blue',
		"MEF_Dox_D2_rep1" = 'purple',
		"MEF_Dox_D2_rep2" = 'purple',
		"MEF_Dox_D7_rep1" = 'red',
		"MEF_Dox_D7_rep2" = 'red',
		"MEF_Dox_D7_Flk1pos_rep1" = 'pink',
		"MEF_Dox_D7_Flk1pos_rep2" = 'pink',
		"EB_Dox_D25_rep1" = 'red',
		"EB_Dox_D25_rep2" = 'red',
		"EB_Dox_D25_Flk1pos_rep1" = 'pink',
		"EB_Dox_D25_Flk1pos_rep2" = 'pink',
		"EB_NoDox_D25_rep1" = 'gray',
		"EB_NoDox_D25_rep2" = 'gray'
	)
	g2bg[group]
}

atacseq_group2pch <- function(group){
	g2pch <- c(
		"MEF_NoDox_rep1" = 21,
		"MEF_NoDox_rep2" = 21,
		"MEF_Dox_D1_rep1" = 21,
		"MEF_Dox_D1_rep2" = 21,
		"MEF_Dox_D2_rep1" = 21,
		"MEF_Dox_D2_rep2" = 21,
		"MEF_Dox_D7_rep1" = 21,
		"MEF_Dox_D7_rep2" = 21,
		"MEF_Dox_D7_Flk1pos_rep1" = 21,
		"MEF_Dox_D7_Flk1pos_rep2" = 21,
		"EB_Dox_D25_rep1" = 3,
		"EB_Dox_D25_rep2" = 3,
		"EB_Dox_D25_Flk1pos_rep1" = 3,
		"EB_Dox_D25_Flk1pos_rep2" = 3,
		"EB_NoDox_D25_rep1" = 3,
		"EB_NoDox_D25_rep2" = 3
	)
	g2pch[group]
}

atacseq_group2col <- function(group){
	g2col <- c(
		"MEF_NoDox_rep1" = 'gray', 
		"MEF_NoDox_rep2" = 'gray',
		"MEF_Dox_D1_rep1" = 'black',
		"MEF_Dox_D1_rep2" = 'black',
		"MEF_Dox_D2_rep1" = 'black',
		"MEF_Dox_D2_rep2" = 'black',
		"MEF_Dox_D7_rep1" = 'black',
		"MEF_Dox_D7_rep2" = 'black',
		"MEF_Dox_D7_Flk1pos_rep1" = 'black',
		"MEF_Dox_D7_Flk1pos_rep2" = 'black',
		"EB_Dox_D25_rep1" = 'red',
		"EB_Dox_D25_rep2" = 'red',
		"EB_Dox_D25_Flk1pos_rep1" = 'pink',
		"EB_Dox_D25_Flk1pos_rep2" = 'pink',
		"EB_NoDox_D25_rep1" = 'green',
		"EB_NoDox_D25_rep2" = 'green'
	)
	g2col[group]
}


# --------------------------------------------------------------------------------
# [2019-05-03] Read Etv2 ChIP-seq peaks in MEF
# --------------------------------------------------------------------------------
read_Etv2_peaks_in_MEF <- function(exclude_promoter = TRUE, exclude_exons = TRUE, width = 2000, minoverlap = 30){

	library(BSgenome.Mmusculus.UCSC.mm10)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)

	peaks_files <- sprintf('analysis/etv2_pioneer/results/%s_summits.bed', c(
		'MEF_Dox_D1_Etv2',
		'MEF_Dox_D2_Etv2',
		'MEF_Dox_D7_Etv2'
	))
	gr0 <- lapply(peaks_files, function(file){
		flog.info(sprintf('reading %s', file))	
		gr <- read.table(file, header = FALSE, sep = '\t')
		gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]))
		gr <- reCenterPeaks(gr, width = width)
		gr <- add.seqinfo(gr, genome = 'mm10')
		gr
	})

	gr <- gr0[[1]]
	for (i in 2:length(gr0)){
		mm <- as.matrix(findOverlaps(gr, gr0[[i]], minoverlap = minoverlap))
		gr <- c(gr, gr0[[i]][!1:length(gr0[[i]]) %in% mm[, 'subjectHits']])
	}

	if (exclude_promoter){
		pmt <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene)
		is_promoter <- 1:length(gr) %in% as.matrix(findOverlaps(gr, pmt))[, 1]
		gr <- gr[!is_promoter]
		flog.info(sprintf('get %d intervals after removing promoter region', length(gr)))
	}

	if (exclude_exons){
		ex <- exons(TxDb.Mmusculus.UCSC.mm10.knownGene)
		is_exon <- 1:length(gr) %in% as.matrix(findOverlaps(gr, ex))[, 1]
		gr <- gr[!is_exon]
		flog.info(sprintf('get %d intervals after removing exon region', length(gr)))
	}
	gr
}


# --------------------------------------------------------------------------------
# [2019-05-03] Read Etv2 ChIP-seq peaks in MEF
# --------------------------------------------------------------------------------
read_Etv2_peaks_in_D1_MEF <- function(width = 200, exclude_promoter = TRUE, exclude_exons = TRUE, MEF_H3K27ac = NULL){

	library(BSgenome.Mmusculus.UCSC.mm10)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)

	peak_file <- sprintf('%s/macs2/MEF_Dox_D1_Etv2_summits.bed', PROJECT_DIR)
	flog.info(sprintf('reading %s', peak_file))	
	gr <- read.table(peak_file, header = FALSE, sep = '\t')
	gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]))
	gr <- reCenterPeaks(gr, width = width)
	gr <- add.seqinfo(gr, genome = 'mm10')

	if (exclude_promoter){
		pmt <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene)
		is_promoter <- 1:length(gr) %in% as.matrix(findOverlaps(gr, pmt))[, 1]
		gr <- gr[!is_promoter]
		flog.info(sprintf('get %d intervals after removing promoter region', length(gr)))
	}

	if (exclude_exons){
		ex <- exons(TxDb.Mmusculus.UCSC.mm10.knownGene)
		is_exon <- 1:length(gr) %in% as.matrix(findOverlaps(gr, ex))[, 1]
		gr <- gr[!is_exon]
		flog.info(sprintf('get %d intervals after removing exon region', length(gr)))
	}

	if (!is.null(MEF_H3K27ac)){
		bw_files <- get_bigwig_files()
		bw_files <- bw_files['MEF_H3K27ac']
		bed_files <- gsub('_treat_pileup.bw', '_peaks.broadPeak', bw_files)
		h3k27ac <- read.table(bed_files, header = FALSE, sep = '\t')
		h3k27ac <- GRanges(seqnames = h3k27ac[, 1], range = IRanges(h3k27ac[, 2], h3k27ac[, 3]))
		h3k27ac <- add.seqinfo(h3k27ac, genome = 'mm10')

		if (MEF_H3K27ac){
			gr <- subsetByOverlaps(gr, h3k27ac)
			flog.info(sprintf('get %d intervals after including H3K27ac ON regions', length(gr)))
		}else{
			gr <- subsetByOverlaps(gr, h3k27ac, invert = TRUE)
			flog.info(sprintf('get %d intervals after excluding H3K27ac ON regions', length(gr)))
		}
	}
	gr

} # read_Etv2_peaks_in_D1_MEF



# --------------------------------------------------------------------------------
# [2019-05-03] Get a list of motifs
# --------------------------------------------------------------------------------
get_motifs <- function(tf = 'Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer(0.967)', exclude_promoter = TRUE, exclude_exons = TRUE){

	library(chromVARmotifs) # https://github.com/GreenleafLab/chromVARmotifs
	library(motifmatchr)
	library(BSgenome.Mmusculus.UCSC.mm10)
	library(TxDb.Mmusculus.UCSC.mm10.knownGene)

	motif.set <- 'homer_pwms'
	data(list = motif.set, package = 'chromVARmotifs')

	# a GR of all genomic regions
	gr <- GRanges(seqnames = seqnames(BSgenome.Mmusculus.UCSC.mm10), range = IRanges(1, seqlengths(BSgenome.Mmusculus.UCSC.mm10)))
	gr <- add.seqinfo(gr, genome = 'mm10')
	se <- SummarizedExperiment(rowRanges = gr)
	motif_ix <- matchMotifs(get(motif.set)[tf], se, genome = 'mm10', out = 'positions')
	gr <- motif_ix[[1]]
	gr <- add.seqinfo(gr, genome = 'mm10')
	flog.info(sprintf('get %d intervals', length(gr)))

	if (exclude_promoter){
		pmt <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene)
		is_promoter <- 1:length(gr) %in% as.matrix(findOverlaps(gr, pmt))[, 1]
		gr <- gr[!is_promoter]
		flog.info(sprintf('get %d intervals after removing promoter region', length(gr)))
	}

	if (exclude_exons){
		ex <- exons(TxDb.Mmusculus.UCSC.mm10.knownGene)
		is_exon <- 1:length(gr) %in% as.matrix(findOverlaps(gr, ex))[, 1]
		gr <- gr[!is_exon]
		flog.info(sprintf('get %d intervals after removing exon region', length(gr)))
	}

	gr
} 

# --------------------------------------------------------------------------------
# [2019-05-03] Shift the ATAC-seq reads
# Modified from https://rdrr.io/bioc/ATACseqQC/src/R/shiftReads.R
# Only isize field is available in mcols(x)
# --------------------------------------------------------------------------------
shiftReads <- function(x, positive=4L, negative=5L){

	strds <- as.character(strand(x)) == "-"
  ns <- ifelse(strds, negative, positive)

	## fix soft-clipping
	cigars <- cigar(x)
#	mcols(x)$seq <- sequenceLayer(mcols(x)$seq, cigars, from="query", to="query-after-soft-clipping")
#	mcols(x)$qual <- sequenceLayer(mcols(x)$qual, cigars, from="query", to="query-after-soft-clipping")
	cigars <- as.character(cigarNarrow(cigars))
	cigars <- cigarQNarrow(cigars, start=ifelse(strds, 1, positive+1), end=ifelse(strds, -negative-1, -1))
	x@cigar <- as.character(cigars)
	x@start <- x@start + attributes(cigars)$rshift
				  
	width.x <- qwidth(x)
	mcols(x)$isize <- sign(mcols(x)$isize) * (abs(mcols(x)$isize) - positive - negative)
	##x end will auto matic change if the cigar changed.
#	mcols(x)$seq[strds] <- DNAStringSet(substr(mcols(x)$seq[strds], 1,width.x[strds]))
#	mcols(x)$qual[strds] <- PhredQuality(substr(mcols(x)$qual[strds], 1,width.x[strds]))
#	mcols(x)$seq[!strds] <- DNAStringSet(substr(mcols(x)$seq[!strds],ns[!strds]+1,nchar(mcols(x)$seq[!strds])))
#	mcols(x)$qual[!strds] <- PhredQuality(substr(mcols(x)$qual[!strds],ns[!strds]+1,nchar(mcols(x)$qual[!strds])))
	x
}


vplot <- function(atac, gr, upstream = 500, downstream = 500){

	wid <- width(gr)[1]

	flog.info(sprintf('there are %d input ATAC-seq read pairs', length(atac)))
	flog.info(sprintf('there are %d input motifs/peaks', length(gr)))

	flog.info(sprintf('resizing motifs to [-%d, +%d]', upstream, downstream))
	gr <- promoters(reCenterPeaks(gr, width = 1), upstream = upstream+floor(wid/2), downstream=downstream+ceiling(wid/2))

	atac <- subsetByOverlaps(atac, gr)  # ATAC-seq peaks overlapping with the motifs
	flog.info(sprintf('there are %d ATAC-seq read pairs overlapping with motifs', length(atac)))

	flog.info(sprintf('adjusting for the cut sites'))
	atac <- GAlignmentPairs(shiftReads(first(atac)), shiftReads(second(atac)))

	x <- as(atac, 'GRanges')
	x$w <- width(x)
	x <- reCenterPeaks(x, width = 1)
	ol <- findOverlaps(x, gr)
	x_query <- x[queryHits(ol)]
	gr_subject <- gr[subjectHits(ol)]
	rel <- getRelationship(x_query, gr_subject)
	rel$FragmentLength <- x_query$w
	rel$distanceToBindingSite <- rel$distanceToStart - upstream - floor(wid/2)
	rel <- rel[order(rel$distanceToStart, rel$FragmentLength), c("distanceToBindingSite", "FragmentLength")]
	rel
}

# https://rdrr.io/bioc/ATACseqQC/src/R/vPlot.R
getRelationship <- function(queryHits, subjectHits){
	if(!inherits(queryHits, "GRanges")) 
		stop("queryHits must be an object of GRanges")
  if(!inherits(subjectHits, "GRanges")) 
		stop("subjectHits must be an object of GRanges")
	strand <- strand(subjectHits)=="-"
	FeatureStart <- as.numeric(ifelse(strand, end(subjectHits), start(subjectHits)))
	FeatureEnd <- as.numeric(ifelse(strand, start(subjectHits), end(subjectHits)))
	PeakStart <- as.numeric(ifelse(strand, end(queryHits), start(queryHits)))
	PeakEnd <- as.numeric(ifelse(strand, start(queryHits), end(queryHits)))
	ss <- PeakStart - FeatureStart
	ee <- PeakEnd - FeatureEnd
	se <- PeakStart - FeatureEnd
	es <- PeakEnd - FeatureStart
	shortestDistance <- apply(cbind(ss, ee, se, es), 1,function(.ele) min(abs(.ele)))
	shortestDistanceToStart <- apply(cbind(ss, es), 1, function(.ele) min(abs(.ele)))
	data.frame(shortestDistance=shortestDistance, ss=ss,distanceToStart=shortestDistanceToStart)
}


read_ATAC_dataset <- function(sample_level = TRUE, touch = TRUE){

	dataset <- 'dataset=Etv2ATAC_version=20190228a'
	d <- read.dataset(dataset, touch = touch)
	d <- transform(d, sra.run.dir = sra.run.dir(run), bam.file = sprintf('%s/%s.dedup.bam', sra.run.result.dir(run), run))
	d <- transform(d, group2 = gsub('_rep\\d$', '', group))

	flog.info('ga_file= *_.rds')
	d <- transform(d, ga_file = gsub('.bam', '.rds', bam.file)) # GAlignments file for each BAM

	flog.info('fragment_size_file = *_.fragment_size.txt')
	d <- transform(d, fragment_size_file = gsub('.bam', '.fragment_size.txt', bam.file))	# global ATAC-seq fragment size file

	for (i in 1:nrow(d)){
		for (f in c('ga_file', 'fragment_size_file')){
			command <- sprintf('touch -h %s', d[i, f])  # adding -h option for symlink (see https://unix.stackexchange.com/questions/63876/changing-the-timestamp-of-a-symlink)
			flog.info(command)
			system(command)
		}
	}

	if (sample_level){
		sp <- split(d[, 'bam.file'], list(d[, 'group2']))
		bam_files <- sprintf('analysis/etv2_pioneer/data/ATAC_%s.bam', names(sp))
		d2 <- data.frame(group = names(sp), bam_file = bam_files)
		rownames(d2) <- d2[, 'group']
		d2 <- transform(d2, ga_file = gsub('.bam', '.rds', bam_file)) # GAlignments file for each BAM
		d2 <- transform(d2, vplot_file = gsub('.bam', '_Etv2peaks_vplot.rds', bam_file))  # Vplot density file for each BAM
		d2 <- transform(d2, vplot_motif_file = gsub('.bam', '_Etv2motif_vplot.rds', bam_file))  # Vplot density file for each BAM
		d2 <- transform(d2, vplot_Etv2_D1_MEF_file = gsub('.bam', '_Etv2peaks_D1_MEF_vplot.rds', bam_file)) # Vplot density file for each BAM
		d2
	}else{
		d
	}

} # read_ATAC_dataset


# --------------------------------------------------------------------------------
# [2019-05-14] Read Etv2 ChIP-seq peaks in MEF, Etv2 peaks from Choi's data
# [2019-05-14] To make it simplier, for Choi's data we only used V5 antibodies.  
# --------------------------------------------------------------------------------
read_Etv2_peaks <- function(width = 500){

	library(BSgenome.Mmusculus.UCSC.mm10)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
	library(ChIPpeakAnno)

	# MEF Etv2 peaks
	peaks_files <- c(
		MEF_Dox_D1 = 'analysis/etv2_pioneer/results/MEF_Dox_D1_Etv2_summits.bed',
		MEF_Dox_D2 = 'analysis/etv2_pioneer/results/MEF_Dox_D2_Etv2_summits.bed',
		MEF_Dox_D7 = 'analysis/etv2_pioneer/results/MEF_Dox_D7_Etv2_summits.bed',
		EB_D4_V5 = 'https://s3.msi.umn.edu/etv2_pioneer/etv2_chipseq_V5_mm10.sorted.bed'
#		EB_D4_pAb = 'https://s3.msi.umn.edu/etv2_pioneer/etv2_chipseq_polyab_mm10.sorted.bed'
	)
	ss <- names(peaks_files)

	gr0 <- lapply(peaks_files, function(file){
		flog.info(sprintf('reading %s', file))	
		gr <- read.table(file, header = FALSE, sep = '\t')
		gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]))
		gr <- reCenterPeaks(gr, width = width)
		gr <- add.seqinfo(gr, genome = 'mm10')
		gr <- unique(gr)
		gr
	})

#	x <- findOverlapsOfPeaks(gr0$EB_D4_V5, gr0$EB_D4_pAb)
#	gr0$EB_D4 <- x$peaklist[['gr0.EB_D4_V5///gr0.EB_D4_pAb']]
#	gr0$EB_D4 <- Reduce('c', x$peaklist)
	gr0$EB_D4 <- gr0$EB_D4_V5

	x <- findOverlapsOfPeaks(gr0$MEF_Dox_D1, gr0$MEF_Dox_D2, gr0$MEF_Dox_D7, gr0$EB_D4)
	print(x$venn_cnt)

	y <- names(x$peaklist)
	G <- x$venn_cnt[, -ncol(x$venn_cnt)]
	colnames(G) <- ss
	G <- G[-1, ] # all zero's

	h <- sapply(x$peaklist,length)
	G <- G[rep(1:nrow(G), h), ]
	gr <- Reduce('c', x$peaklist)
	mcols(gr)$group <- G
	class(G) <- 'numeric'
	mcols(gr)$cluster <- Reduce('paste0', lapply(1:ncol(G), function(i) G[, i]))
	metadata(gr)$ol <- x
	gr

} # read_Etv2_peaks 


get_bigwig_color <- function(){
	cols <- c(
		'MEF_Dox_D1_Etv2' = 'white',
		'MEF_Dox_D2_Etv2' = 'white',
		'MEF_Dox_D7_Etv2' = 'white',
		'EB_D4_Etv2_pAb' = 'white',
		'EB_D4_Etv2_V5' = 'white',
		'MEF_nucleosome' = 'green',
		'MEF_Brg1' = 'red',
		'MEF_H3' = 'yellow',
		'MEF_H3K9me3' = 'lightblue',
		'MEF_H3K27me3' = 'olivedrab1',

		'MEF_ATAC' = 'skyblue',# public ATAC-seq data in MEF

		'MEF_NoDox_ATAC' = 'skyblue',
		'MEF_D1_ATAC' = 'skyblue',
		'MEF_D2_ATAC' = 'skyblue',
		'MEF_D7_ATAC' = 'skyblue',
		'MEF_D7_Flk1pos_ATAC' = 'skyblue',

		'EB_NoDox_D25' = 'skyblue',
		'EB_Dox_D25' = 'skyblue',
		'EB_Dox_D25_Flk1pos' = 'skyblue',

		'MEF_H3K36me3' = 'brown',
		'MEF_H3K9ac' = 'hotpink',
		'MEF_H3K79me2' = 'brown',
		'MEF_H3K4me2' = 'mediumpurple',
		'MEF_H3K4me1' = 'lightgoldenrod',
		'MEF_Hdac1' = 'orchid',
		'MEF_H3.3' = 'rosybrown1',
		'MEF_P300' = 'gold',

		'MEF_H3K27ac' = 'yellow',
		'MEF_NoDox_H3K27ac' = 'yellow',
		'MEF_D1_H3K27ac' = 'yellow',
		'MEF_D2_H3K27ac' = 'yellow',
		'MEF_D7_H3K27ac' = 'yellow',

		# in house repeat level ATAC-seq data
		'MEF_NoDox_rep1' = 'skyblue',
		'MEF_NoDox_rep2' = 'skyblue',
		'MEF_Dox_D1_rep1' = 'skyblue',
		'MEF_Dox_D1_rep2' = 'skyblue',
		'MEF_Dox_D2_rep1' = 'skyblue',
		'MEF_Dox_D2_rep2' = 'skyblue',
		'MEF_Dox_D7_rep1' = 'skyblue',
		'MEF_Dox_D7_rep2' = 'skyblue',
		'MEF_Dox_D7_Flk1pos_rep1' = 'skyblue',
		'MEF_Dox_D7_Flk1pos_rep2' = 'skyblue',
		'EB_NoDox_D25_rep1' = 'skyblue',
		'EB_NoDox_D25_rep2' = 'skyblue',
		'EB_Dox_D25_rep1' = 'skyblue',
		'EB_Dox_D25_rep2' = 'skyblue',
		'EB_Dox_D25_Flk1pos_rep1' = 'skyblue',
		'EB_Dox_D25_Flk1pos_rep2' = 'skyblue'
	)
	cols
} # get_bigwig_color


read_H3K27ac_dataset <- function(sample_level = FALSE, touch = TRUE){

	dataset <- 'dataset=Etv2ChIPseq_version=20190307a'
	d <- read.dataset(dataset, touch = touch)
	d <- transform(d, sra.run.dir = sra.run.dir(run), bam.file = sprintf('%s/%s.dedup.bam', sra.run.result.dir(run), run))
	d <- transform(d, group2 = gsub('_rep\\d$', '', group))
	d <- d[d$seqtype %in% c('H3K27ac', 'input'), ]

	dataset <- 'dataset=PioneerEB_version=20190501a'
	d2 <- read.dataset(dataset, touch = touch)
	d2 <- transform(d2, sra.run.dir = sra.run.dir(run), bam.file = sprintf('%s/%s.dedup.bam', sra.run.result.dir(run), run))
	d2 <- transform(d2, group2 = gsub('_rep\\d$', '', group))

	rbind(d, d2)

} # read_H3K27ac_dataset


# --------------------------------------------------------------------------------
# [2019-05-16] Read ATAC-seq peaks of our ATAC-seq data from MEF and ES/EB
# --------------------------------------------------------------------------------
read_ATAC_peaks <- function(width = 200, log10qvalue = 0, exclude_promoter = TRUE, exclude_exons = TRUE){

	library(BSgenome.Mmusculus.UCSC.mm10)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
	library(ChIPpeakAnno)

	peaks_files <- c(
#		'MEF_ATAC' = 'analysis/etv2_pioneer/data/Chronis/ATAC_summits.bed',
#		'MEF_NoDox_ATAC' = 'ATAC_MEF_NoDox_summits.bed',
		'MEF_D1_ATAC' = 'ATAC_MEF_Dox_D1_summits.bed',
#		'MEF_D2_ATAC' = 'ATAC_MEF_Dox_D2_summits.bed',
#		'MEF_D7_ATAC' = 'ATAC_MEF_Dox_D7_summits.bed',
		'MEF_D7_Flk1pos_ATAC' = 'ATAC_MEF_Dox_D7_Flk1pos_summits.bed',
#		'EB_NoDox_D25' = 'ATAC_EB_NoDox_D25_summits.bed',
#		'EB_Dox_D25' = 'ATAC_EB_Dox_D25_summits.bed'
		'EB_Dox_D25_Flk1pos' = 'ATAC_EB_Dox_D25_Flk1pos_summits.bed'
	)
	peaks_files <- sprintf('%s/macs2/%s', PROJECT_DIR, peaks_files)

	ss <- names(peaks_files)

	gr0 <- lapply(peaks_files, function(file){
		flog.info(sprintf('reading %s', file))	
		gr <- read.table(file, header = FALSE, sep = '\t')
		gr <- GRanges(seqnames = gr[, 1], range = IRanges(gr[, 2], gr[, 3]), log10qvalue = gr[, 5])

		flog.info(sprintf('use the peaks with log10 of qvalue > %.3f', log10qvalue))
		gr <- gr[mcols(gr)$log10qvalue > log10qvalue]	# get rid of unclear peaks

		gr <- reCenterPeaks(gr, width = width)
		gr <- add.seqinfo(gr, genome = 'mm10')
		gr <- unique(gr)
		gr
	})

	gr <- gr0[[1]]

	for (i in 2:length(gr0)){
		x <- findOverlapsOfPeaks(gr, gr0[[i]])
		gr <- Reduce('c', x$peaklist)
	}

	G <- do.call('cbind', lapply(1:length(gr0), function(i) 1:length(gr) %in% as.matrix(findOverlaps(gr, gr0[[i]]))[, 1]))
	colnames(G)  <- ss
	mcols(gr)$group <- G
	class(G) <- 'numeric'
	mcols(gr)$cluster <- Reduce('paste0', lapply(1:ncol(G), function(i) G[, i]))

	if (exclude_promoter){
		pmt <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene)
		is_promoter <- 1:length(gr) %in% as.matrix(findOverlaps(gr, pmt))[, 1]
		gr <- gr[!is_promoter]
		flog.info(sprintf('get %d intervals after removing promoter region', length(gr)))
	}

	if (exclude_exons){
		ex <- exons(TxDb.Mmusculus.UCSC.mm10.knownGene)
		is_exon <- 1:length(gr) %in% as.matrix(findOverlaps(gr, ex))[, 1]
		gr <- gr[!is_exon]
		flog.info(sprintf('get %d intervals after removing exon region', length(gr)))
	}
	gr

} # read_ATAC_peaks


# --------------------------------------------------------------------------------
# [2019-05-17] Get the PWMs which genes are expressed in either D1 or MEF according to the scRNA-seq
# --------------------------------------------------------------------------------
get_PWM_in_MEF <- function(){

	library(motifmatchr)
	library(BSgenome.Mmusculus.UCSC.mm10)
	library(biomaRt)

	motif.set <- 'mouse_pwms_v2'
	data(list = motif.set, package = 'chromVARmotifs')
	pwms <- get(motif.set)

	d <- data.frame(id = 1:length(pwms), gene_id = names(pwms))
	d <- transform(d, is_ensembl_gene_id = grepl('ENSMUSG', gene_id), is_refseq_protein = grepl('NP_\\d+', gene_id) | grepl('XP_\\d+', gene_id))
	d$ensembl_gene_id <- rep(NA, nrow(d))
	d$ensembl_gene_id[d$is_ensembl_gene_id] <- gsub('(ENSMUSG\\d+).+', '\\1', d$gene_id[d$is_ensembl_gene_id])
	d <- d[d$is_ensembl_gene_id, ]

	flog.info('reading scRNA-seq in MEF')
	dataset <- 'dataset=Etv2scRNAseq_version=20190416b'
	se <- readRDS(sprintf('analysis/datasets/%s.rds', dataset))

	flog.info(sprintf('reading %s', dataset))
	m <- colData(se)$group %in% c('MEF_NoDox', 'MEF_Dox_D1')
	X <- assays(se)$counts[, m]
	gs <- rowData(se)$id[rowSums(X > 1) >= 50]	# genes expressed in Dox D1 or MEF

	d <- d[d$ensembl_gene_id %in% gs, ]
	mart <- useMart(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl', host = 'www.ensembl.org')
	bm <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', values = d[, 'ensembl_gene_id'], mart = mart)
	d <- merge(d[, c('id', 'ensembl_gene_id')], bm, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')
	pwms <- pwms[d[, 'id']]
	names(pwms) <- d$external_gene_name
	pwms

} # get_PWM_in_MEF


# --------------------------------------------------------------------------------
# [2019-05-21] Get PWMs from mouse_pwms_v2
# Map all the gene IDs to standard gene symbol
# --------------------------------------------------------------------------------
get_mouse_pwms_v2 <- function(){

	library(motifmatchr)
	library(BSgenome.Mmusculus.UCSC.mm10)
	library(biomaRt)

	motif.set <- 'mouse_pwms_v2'
	data(list = motif.set, package = 'chromVARmotifs')
	pwms <- get(motif.set)

	d <- data.frame(id = 1:length(pwms), gene_id = names(pwms))
	d <- transform(d, is_ensembl_gene_id = grepl('ENSMUSG', gene_id), is_refseq_protein = grepl('NP_\\d+', gene_id) | grepl('XP_\\d+', gene_id))
	d$ensembl_gene_id <- rep(NA, nrow(d))
	d$ensembl_gene_id[d$is_ensembl_gene_id] <- gsub('(ENSMUSG\\d+).+', '\\1', d$gene_id[d$is_ensembl_gene_id])
	d <- d[d$is_ensembl_gene_id, ]
	mart <- useMart(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl', host = 'www.ensembl.org')
	bm <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', values = d[, 'ensembl_gene_id'], mart = mart)
	d <- merge(d[, c('id', 'ensembl_gene_id')], bm, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')
	pwms <- pwms[d[, 'id']]
	names(pwms) <- d$external_gene_name
	pwms

} # get_PWMs




# ------------------------------------------------------------------------------------------------------------------
# [2019-05-21] Generate PFM input for PIQ; the PFM must be stored as the JASPAR format according to 
# https://bitbucket.org/thashim/piq-single/src/master/
# Note: the mouse_pwms_v2 dataset from chromVARmotifs package is represented as the PWM matrix while the PIQ requires the PFM input
# There appears no ready-to-use function to convert PWM to PFM. 
# ------------------------------------------------------------------------------------------------------------------
write_jaspar_for_PIQ <- function(pwms, file, n = 100){

	library(BSgenome.Mmusculus.UCSC.mm10)
	library(universalmotif)	# write_jaspar 
	library(TFBSTools)	# PFMatrix

	gr <- GRanges(seqnames = seqnames(BSgenome.Mmusculus.UCSC.mm10), range = IRanges(1, seqlengths(BSgenome.Mmusculus.UCSC.mm10)))
	gr <- gr[seqnames(gr) == 'chr1']	# only use chromsome 1
	flog.info(sprintf('running matchMotifs of %d motifs on chromosome 1', length(pwms)))
	motif_ix <- matchMotifs(pwms, gr, genome = 'mm10', out = 'positions')
	y <- bplapply(motif_ix, function(x) x[strand(x) == '+'])
	y <- bplapply(y, function(x) x[sample.int(length(x), n)])
	flog.info(sprintf('get the DNA sequences'))
	S <- bplapply(y, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))
	S <- bplapply(S, function(x) do.call('rbind', strsplit(as.character(x), '')))
	M <- bplapply(S, function(x) apply(x, 2, function(z) table(factor(z, c('A', 'C', 'G', 'T')))))
	pfm_list <- bplapply(1:length(pwms), function(i) PFMatrix(ID = names(pwms)[i], name = names(pwms)[i], matrixClass = 'Unknown', strand= '+', bg = pwms[[i]]@bg, tags = list(family = 'Unknown'), profileMatrix = M[[i]]))
	pfm_list <- lapply(pfm_list, function(x) convert_motifs(x))
	flog.info(sprintf('writing %s', file))
	write_jaspar(pfm_list, file)

} #  write_jaspar_for_PIQ

matrix2normalizedMatrix <- function(X, center = 300, extend = 250){

	X <- (X[, 301:600] + X[, 300:1]) / 2
	X <- as.matrix(Diagonal(x = 1 / rowSums(X)) %*% X)
	X <- t(apply(X, 1, default_smooth_fun))
	X <- X[, 1:extend]
	attr(X, 'upstream_index') <- 1
	attr(X, 'target_index') <- 1:72
	attr(X, 'downstream_index') <- 73:250
	attr(X, 'extend') <- c(0, 250)
	attr(X, 'target_is_single_point') <- FALSE
	class(X) <- c("normalizedMatrix", "matrix")
	X

}


