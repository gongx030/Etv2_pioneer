# ---------------------------------------------------------------------
# [2019-11-18] chromVAR analysis of TFs dynamics of Etv2 reprogramming ATAC-seq data
# ---------------------------------------------------------------------
bed_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/ATAC_peaks_Etv2_reprogramming.bed'
devtools::load_all('packages/compbio'); peaks <- macs2.read_summits(bed_file)


bw_files <- c(
	'MEF_NoDox' = '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_NoDox_treat_pileup.bw',
	'MEF_Dox_D1' = '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D1_treat_pileup.bw',
	'MEF_Dox_D2' = '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D2_treat_pileup.bw',
	'MEF_Dox_D7' = '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D7_treat_pileup.bw',
	'MEF_Dox_D7_Flk1pos' = '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/MEF_Dox_D7_Flk1pos_treat_pileup.bw'
)
X <- do.call('cbind', mclapply(bw_files, function(bw_file){
	flog.info(sprintf('reading %s', bw_file))
	ga <- rtracklayer::import(bw_file, format = 'BigWig', which = reduce(peaks))
	ga <- add.seqinfo(ga, 'mm10')
	cvg <- coverage(ga, weight = as.numeric(mcols(ga)$score))
	sum(cvg[peaks])
}, mc.cores = 2))


library(chromVARmotifs) # https://github.com/GreenleafLab/chromVARmotifs
library(chromVAR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)
library(BiocParallel)
library(SummarizedExperiment)
register(MulticoreParam(4)) # Use 8 cores
motif.set <- 'mouse_pwms_v2'
data(list = motif.set, package = 'chromVARmotifs')
se <- SummarizedExperiment(assays = SimpleList(counts = X), rowRanges = peaks)
se <- addGCBias(se, genome = BSgenome.Mmusculus.UCSC.mm10)
motif_ix <- matchMotifs(get(motif.set), peaks, genome = 'mm10')
dev <- computeDeviations(object = se, annotations = motif_ix)
v <- computeVariability(dev)


