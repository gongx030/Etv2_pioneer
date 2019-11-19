
# ----------------------------------------------------------------------------
# [2019-04-12] Comparing two Brg1 ChIP-seq peaks and see the consistency between them
# If they are largely consistent, this will reduce the difficulities of re-doing the Brg1 ChIP-seq
# [2019-09-06] Update the analysis using latest Brg1 ChIP-seq data
# ----------------------------------------------------------------------------
devtools::load_all('packages/compbio')
bed_files <- c(
	'/panfs/roc/scratch/gongx030/datasets/dataset=Chronis_version=20190405a/Brg1_summits.bed',
	'/panfs/roc/scratch/gongx030/datasets/dataset=Alver_version=20190407a/Brg1_summits.bed',
	'/panfs/roc/scratch/gongx030/datasets/dataset=Brg1_version=20190820a/MEF_NoDox_Brg1_summits.bed',
	'/panfs/roc/scratch/gongx030/datasets/dataset=Brg1_version=20190820a/MEF_Dox_D1_Brg1_summits.bed'
)
grs <- lapply(bed_files, function(bed_file) macs2.read_summits(bed_file, score_threshold = -log10(0.001)))
grs <- lapply(grs, function(gr) resize(gr, fix = 'center', width = 500))
library(ChIPpeakAnno); ol <- findOverlapsOfPeaks(grs[[1]], grs[[2]], grs[[3]], grs[[4]])
makeVennDiagram(ol)


# ----------------------------------------------------------------------------
# [2019-11-18] Preprocessing the Etv2 scRNA-seq data
# Remove the potential doublet using scrublet
# ----------------------------------------------------------------------------
devtools::load_all('packages/compbio')
dataset <- 'dataset=Etv2scRNAseq_version=20190416a'
se <- readRDS(sprintf('%s/filtered.rds', dataset_dir(dataset)))
i <- colSums(assays(se)$counts) > 5000
se <- se[, i]
saveRDS(se, sprintf('%s/filtered_5000.rds', dataset_dir(dataset)))
i <- rowSums(assays(se)$counts > 0) >= 5
se <- se[i, ]
saveRDS(se, sprintf('%s/filtered_5000_fivegenes.rds', dataset_dir(dataset)))
devtools::load_all('packages/compbio'); se <- scrublet(se, colData(se)$group, expected_doublet_rate = 0.06, min_counts = 2, min_cells = 3, min_gene_variability_pctl = 85, n_prin_comps = 30L)
se <- se[, !colData(se)$predicted_doublets]
saveRDS(se, sprintf('%s/filtered_5000_fivegenes_noDoublets.rds', dataset_dir(dataset)))


# ----------------------------------------------------------------------------
# [2019-04-16] Preprocessing the Etv2 scRNA-seq data
# Need to run on lab queue, will fail on the k40
# [2019-11-18] Preprocessing scRNA-seq using Seurat (3.1.1)
# ----------------------------------------------------------------------------
se_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2scRNAseq_version=20190416a/filtered_5000_fivegenes_noDoublets.rds'
se <- readRDS(se_file)
library(BiocParallel)
register(MulticoreParam(16))
devtools::load_all('packages/compbio'); se2 <- scran_preprocess(se)
se_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2scRNAseq_version=20190416a/filtered_5000_fivegenes_noDoublets_scran.rds'
saveRDS(se2, se_file)


# ----------------------------------------------------------------------------
# [2019-04-17] scVI analysis of the scRNA-seq of Etv2 induction in mouse
# Run this on k40
# [2019-04-23] SCRAN gives less HVGs
# [2019-11-18] Updated analysis
# ----------------------------------------------------------------------------
library(SummarizedExperiment)
se_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2scRNAseq_version=20190416a/filtered_5000_fivegenes_noDoublets_scran.rds'
se <- readRDS(se_file)
n <- rowData(se)$FDR < 0.05
latent <- 10
devtools::load_all('packages/compbio'); se2 <- scvi(se[n, ], n_latent = latent)
colData(se)$latent <- colData(se2)$latent
se_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2scRNAseq_version=20190416a/filtered_5000_fivegenes_noDoublets_scran_scVI.rds'
saveRDS(se, se_file)

# ----------------------------------------------------------------------------
# [2019-04-23] TSNE of the Etv2 scRNA-seq on MEF
# TSNE runs on the latent space obtained by scVI
# [2019-11-18] Updated visualization using umap
# ----------------------------------------------------------------------------
library(SummarizedExperiment)
library(SingleCellExperiment)
se_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2scRNAseq_version=20190416a/filtered_5000_fivegenes_noDoublets_scran_scVI.rds'
se <- readRDS(se_file)

library(umap); set.seed(1); y <- umap(colData(se)$latent)$layout
colData(se)$umap <- y
#library(Rtsne); set.seed(1); y_tsne <- Rtsne(colData(se)$latent)$Y
group2bg <- c(
	'MEF_Dox_D1' = 'black', 
	'MEF_NoDox' = 'blue', 
	'MEF_Dox_D2' = 'purple', 
	'MEF_Dox_D7a' = 'red', 
	'MEF_Dox_D7b' = 'pink'
)
dev.new(width = 5, height = 5); par(mar = c(2, 2, 2, 2))
bg <- group2bg[colData(se)$group]
plot(colData(se)$umap, cex = 0.5, pch = 21, bg  = bg, col = bg, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')

se_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2scRNAseq_version=20190416a/filtered_5000_fivegenes_noDoublets_scran_scVI_umap.rds'
saveRDS(se, se_file)

# ----------------------------------------------------------------------------
# [2019-04-28] KNN graph and Louvain clustering
# [2019-11-18] Updated analysis based on umap results
# ----------------------------------------------------------------------------
library(igraph)
library(FNN)
devtools::load_all('packages/compbio')
se_file <- '/panfs/roc/scratch/gongx030/datasets/dataset=Etv2scRNAseq_version=20190416a/filtered_5000_fivegenes_noDoublets_scran_scVI_umap.rds'
se <- readRDS(se_file)
k <- 500
knn <- get.knn(colData(se)$umap, k = k)
knn <- data.frame(from = rep(1:nrow(knn$nn.index), k), to = as.vector(knn$nn.index), weight = 1/(1 + as.vector(knn$nn.dist)))
g <- graph_from_data_frame(knn, directed = FALSE)
g <- simplify(g)
lc <- cluster_louvain(g)
clust <- as.numeric(as.factor(membership(lc)))
G_max <- max(unique(clust))
colData(se)$cluster <- clust
dev.new(width = 5, height = 5); par(mar = c(2, 2, 2, 2))
library(RColorBrewer); plot(colData(se)$umap, col = colorRampPalette(brewer.pal(11,'Spectral'))(G_max)[clust], pch = 16, asp = 1, xaxt = 'n', yaxt = 'n', main = 'cluster', xlab = '', ylab = '')
y_centers <- do.call('rbind', lapply(1:G_max, function(i) apply(colData(se)$umap[clust == i, ], 2, median)))
text(y_centers[, 1], y_centers[, 2], 1:G_max, cex = 1.25)



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
library(mclust); set.seed(5); mc <- Mclust(colData(se)$umap, G = 2:G_max)

colData(se)$GMM <- mc$classification
library(RColorBrewer); plot(y_tsne, col = colorRampPalette(brewer.pal(11,'Spectral'))(G_max)[mc$classification], pch = 16, asp = 1, xaxt = 'n', yaxt = 'n', main = 'cluster', xlab = '', ylab = '')
set.seed(1); y_centers <- do.call('rbind', lapply(1:G_max, function(i) y_tsne[sample(which(mc$classification == i), 1), ]))
text(y_centers[, 1], y_centers[, 2], 1:G_max, cex = 3)
table(colData(se)$group, mc$classification)



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


