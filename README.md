### The code and datasets for the manuscript: 

# ETV2 functions as a pioneer factor to regulate and reprogram the endothelial lineage

## Processed datasets

| | Dataset | Format | Script | 
| --- | --- | --- | --- | 
| ATAC-seq counts | [link](https://s3.msi.umn.edu/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/all_ATAC.rds) | SummarizedExperiment | [R](ATAC_seq_preprocess.Rmd) |
| scRNA-seq of Etv2 induced reprogramming | [link](https://s3.msi.umn.edu/gongx030/etv2_pioneer/data/processed_Etv2_scRNAseq.rds) | SummarizedExperiment | [R](scRNA_seq_preprocess.Rmd) |
| Bulk RNA-seq of EB differentiation induced by Etv2 Dox | [link](https://s3.msi.umn.edu/gongx030/datasets/dataset=Etv2RNA-seq_version=20190909a/se.rds) | SummarizedExperiment | |
| A union set of Etv2 ChIP-seq peaks in ES/EB and MEFs | [link](https://s3.msi.umn.edu/gongx030/datasets/dataset=Etv2PioneerChIPseq_version=20191203a/all_Etv2_peaks.rds) | SummarizedExperiment | [R](generate_union_Etv2_peakset.ipynb) |
| Preprocessing Etv2 ChIP-seq, H3K27ac ChIP-seq and Brg1 ChIP-seq in ES/EB and MEF | [link](https://docs.google.com/spreadsheets/d/1UWiduM3Pv-GsVGmfxFApnyVBI1THMR8n8wHg5st3b5c/edit?usp=sharing) | | [R](ChIP_seq_preprocess.Rmd) |
| Bash script for preprocessing the Brg1 KO H3K27ac ChIP-seq in MEF | | | [bash](ChIP-seq_Brg1_KO_H3K27ac_preprocess.sh) |
| Fastq datasets for ChIP-seq, ATAC-seq, bulkRNA-seq, scRNA-seq and NOMeseq | [link](https://docs.google.com/spreadsheets/d/1T3c2zDoSaZWDznwXzm1ZYgpRdAG76lofFHX2Y6HTEAY/edit?usp=sharing) | | |

## Notebooks

|  | Figures | Colab link | 
| --- | --- | --- | 
| Identifying Candidates for shRNA Knock Down testing of additional SWI/SNF Factors  | Extended Data Fig. 9 | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/Etv2_Project_SWI_SNF_shRNA_Candidate_gene_Expression_Testing.ipynb) | 
| Differential Expressing testing candidate genes for knockdown in MEFs | Extended Data Fig. 12a | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/Etv2_Project_Candidate_gene_Expression_Testing.ipynb) |
| scRNA-seq of Etv2 induced reprogramming | Figure 1c-1f <br> Extended Data Fig.  2 <br> Extended Data Fig. 4 | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/scRNA_seq.ipynb) | 
| ATAC-seq and combined ATAC-seq/RNA-seq of <br> Etv2 induced MEF reprogramming and ES/EB differentiation | Figure 1h <br> Figure 1i <br> Extended Data Fig.  6 | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/ATAC_analysis.ipynb) | 
| Overlap of Etv2 ChIP-seq peaks in MEF and EB | Extended Data Fig.  8a and 8b | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_ChIP_seq_peaks.ipynb) | 
| Pathway analysis of early, late and sustained Etv2 peaks | Extended Data Fig.  8c | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Pathway_Etv2_peaks.ipynb) |
| Etv2 motifs in early, late and sustained Etv2 peaks | Figure 3a and 3b | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_motifs_in_early_Etv2_peaks.ipynb) |
| Early, late and sustained Etv2 peaks in MEF | Figure 3g | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/early_Etv2_peaks_in_MEF.ipynb) |
| Early, late and sustained Etv2 peaks in EB | Figure 3h | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/early_late_sustained_Etv2_peaks_in_EB.ipynb) |
| Heatmap of ChIP-seq of sustained ETV2 binding sites in EB differentiation, and ATAC-seq of Control and Brg1 KO EBs at day 4 | Figure 5e | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Enriched_heatmap_of_Etv2_chip_seq_data_and_Brg1_floxed_Brg1_KO.ipynb) |
| D7 Brg1 KD ATAC-seq | Figure 4c | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1_KD_sustained_Etv2_peaks.ipynb) |
| chromVAR analysis of Brg1 KD ATAC-seq at D7| Figure 5a and 5b | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/chromVAR_Brg1_KD_ATAC_D7.ipynb) |
| D7 scRNA-seq including Brg1 KD samples | Figure 4c-4f | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1KD_scRNA_seq_D7.ipynb) | 
| D0 scRNA-seq including Brg1 KD samples | Extended Data Fig.  13a-13g | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1KD_scRNA_seq_D0.ipynb) |
| D0 ChIP-seq H3K27ac in Brg1 KD samples | Extended Data Fig.  13h | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/H3K27ac_Chip_seq_Analysis.ipynb) |
| Lmo2 binding sites | Figure 2g | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Lmo2_track.ipynb) |
| Enriched TFBS in commonly <br> up- and down-regulated genes <br> in EB and MEFs on Etv2 induction | Extended Data Fig.  5i, 5j, 5k | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/TFBS_in_commonly_regulated_genes.ipynb) |
| NFR vs NOR at EB D2.5 | Figure 2d and 2b <br> Extended Data Fig.  7b-7d | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_ChIP_seq_analysis.ipynb) |  |  |
| NFR vs NOR at Etv2 binding sites at D1 MEF reprogramming | Figure 2c | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_binding_D1_MEF.ipynb) |
| NFR vs NOR at Etv2 binding sites at D1 MEF reprogramming <br> (including additional epigenetic marks) | Extended Data Fig.  6a | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_binding_D1_MEF_extended.ipynb) |
| D7 Brg1 KD ATAC-seq (NOR vs NFR) | Figure 4a | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1_KD_sustained_Etv2_peaks_NOR_NFR.ipynb) |
| We look at the enrichment of signals +/-1kb of Etv2 binding sites using Etv2 ChIP-seq 3h and 12h and Brg1 Control and KO ATAC-seq at Day 4. | Extended Data Fig.  15h | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Enriched_heatmap_of_Etv2_chip_seq_data_and_Brg1_floxed_Brg1_KO.ipynb) |
| chromVAR analysis for the change in chromatin accessibililty of Transcription factors at Brg1 Control and KO ATAC-seq at Day 4. | Extended Data Fig.  15g | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/chromVAR_analysis_of_Brg1_floxed_D4_and_Brg1_KO.ipynb) |
| NOMeseq analysis D0 and D1, MNase and Etv2 D1 enriched heatmap showinga few regions of the NOR cluster becoming NFR | Figure 2e| [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/2e_NOME_Enriched_Heatmap_with_MNase_and_Etv2_v1.ipynb) |
| Enriched heatmap showing enrihment at early, late and sustained peaks for MEF Etv2 ChIP and Brg1 KD ChIP data | Extended Data Fig.  14b | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1_KD_ChIP_data_Enriched_Heatmap_at_Etv2_Binding_Sites.ipynb) |
| V-plots of Elk3 motif centric regions in D4 WT EB and D4 Brg1 KO EB | Figure 5f | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Elk3_vplots.ipynb) |
| Motif analysis of Etv2 ChIP-seq peaks | Figure 2f <br> Extended Data Fig. 7e | [R](find_de_novo_motifs_Etv2_chipseq_peaks.Rmd) [R](diff_Etv2_motifs_between_NFR_and_nucleosome.Rmd) |
| Dynamics of Etv2 peak centric V-plots during reprogramming | Extended Data Fig. 10d-10j | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Nucleosome_phasing.ipynb) |
| Fragment size distribution between Flk1+ cells and unsorted mixture cells | Extended Data Fig. 10a-10c | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/FLk1plus_mix.ipynb) |
| Upregulated and downregulated overlapping in MEFs and EBs | Supplemenatary table 4 <br> Supplementary table 5 | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/ATAC_analysis_list_of_upregulated_and_downregulated_genes.ipynb) |

## Main Figures

### Figure 1

|  | Figures | Colab link | 
| --- | --- | --- | 
| The UMAP plot for the scRNA-seq of 948 undifferentiated MEFs, 3,539 reprogrammed cells at 24 hrs, 2,936 cells at 48 hrs and 7,202 cells at 7 days and 827 FLK1+/KDR cells at 7 days post-induction of ETV2 in MEFs | Fig. 1c | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/scRNA_seq.ipynb) | 
| The UMAP plot showing cell clusters from k-means clustering that identified seven distinct cell clusters | Fig. 1d | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/scRNA_seq.ipynb) | 
| The expression profiles of ETV2 and FLK1/KDR | Fig. 1e | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/scRNA_seq.ipynb) | 
| The volcano plot of genes differentially expressed between cluster 1 and cluster 7 | Fig. 1f | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/scRNA_seq.ipynb) |
| The PCA of the variations of transcription factor (TF) associated chromatin accessibility of the ATAC-seq of MEF reprogramming (MEFs, 24 hrs, 48 hrs and 7 days post-induction) and EB differentiation (2.5 days and 3 hrs post induction) | Fig. 1h | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/ATAC_analysis.ipynb) |
| The 31 TF expression levels and motif associated chromatin accessibility consistently showed directional change in both EBs and MEFs (13 up-regulated TFs and 18 down-regulated TFs) | Fig. 1h | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/ATAC_analysis.ipynb) |

### Figure 2

|  | Figures | Colab link | 
| --- | --- | --- | 
|  | Fig. 2a | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_ChIP_seq_analysis.ipynb) | 
| The genomic distribution of EB specific, MEF specific and common ETV2 peaks. The EB and MEF specific ETV2 peaks were more likely distributed at the distal intergenic regions | Fig. 2b | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_ChIP_seq_analysis.ipynb) | 
| The heatmap shows the read density of MNase-seq, BRG1 ChIP-seq and H3K27ac ChIP-seq in MEFs, surrounding 131,001 ETV2 ChIP-seq peaks at 24 hrs post-induction during MEF reprogramming. The ETV2 peaks were divided into four quartiles based on the mean MNase-seq signals of the central 200-bp region | Fig. 2c | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_binding_D1_MEF.ipynb) |
| The heatmap shows the ratio of NFR / nucleosome read density, read density of BRG1 ChIP-seq and H3K27ac ChIP-seq in EBs (day 2.5), surrounding 18,024 ETV2 ChIP-seq at 3 hrs post-induction. The ETV2 peaks were divided into NFR (5,291 peaks) and nucleosome (8,843 peaks) groups according to the local V-plot and fragment size profiles of ATAC-seq day 2.5 EBs without ETV2 induction | Fig. 2d | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_ChIP_seq_analysis.ipynb) |
| The heatmap generated using NOMe-seq shows among 5,320 ETV2 binding sites that were nucleosome occupied at undifferentiated MEF, 4,744 (89.1%) became significantly nucleosome-free (NFR) while 576 (10.9%) stayed NOR at D1 of reprogramming | Fig. 2e | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/2e_NOME_Enriched_Heatmap_with_MNase_and_Etv2_v1.ipynb) | 
| Sequence motif analysis by DREME and CentriMo identified a common GGAAAT motif that were significantly more enriched in NFR regions compared with the nucleosomes in both MEFs and EBs  | Fig. 2f | [R](find_de_novo_motifs_Etv2_chipseq_peaks.Rmd)  |
| A region upstream of Lmo2 that was highly enriched for nucleosomes in both cell types, as measured by MNase-seq and ATAC-seq was selected to perform in vitro nucleosomal binding assays  | Fig. 2g | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Lmo2_track.ipynb)  |

### Figure 3

|  | Figures | Colab link | 
| --- | --- | --- | 
| The average partial Etv2 motif scores in upstream/downstream 500 bp regions surrounding the binding summit of the "early', "late" and "sustained" Etv2 peaks in EB differentiation and MEF reprogramming | Fig. 3a | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_motifs_in_early_Etv2_peaks.ipynb) | 
| The percent of the "early', "late" and "sustained" Etv2 peaks in EBs and MEFs include partial Etv2 motifs in upstream/downstream 50bp regions surrounding the binding summits | Fig. 3b | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_motifs_in_early_Etv2_peaks.ipynb) | 
| The heatmap shows the fold enrichment of ETV2 ChIP-seq, BRG1 ChIP-seq and H3K27ac ChIP-seq at 24 hrs, 48 hrs and 7 days post-induction of ETV2 | Fig. 3g | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/early_Etv2_peaks_in_MEF.ipynb) | 
| The heatmap shows the fold enrichment of ETV2 ChIP-seq, BRG1 ChIP-seq and H3K27ac ChIP-seq, at 3 hrs and 12 hrs post-induction of ETV2 in day 2 EBs | Fig. 3h | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/early_late_sustained_Etv2_peaks_in_EB.ipynb) | 

### Figure 4

|  | Figures | Colab link | 
| --- | --- | --- | 
|  VAE to 18,214 Etv2 ChIP-seq peaks during MEF reprogramming and identified six clusters of V-plot according to the central fragment size distribution | Fig. 4a | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Nucleosome_phasing.ipynb) |
|  The six clusters included three types of V-plots where the central Etv2 sites were nucleosome free (C1, C3 and C4), and three types of V-plots where the central Etv2 sites were nucleosome occupied (C2, C5, and C6), represented by aggregated V-plot | Fig. 4b | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Nucleosome_phasing.ipynb) |
|  The six clusters included three types of V-plots where the central Etv2 sites were nucleosome free (C1, C3 and C4), and three types of V-plots where the central Etv2 sites were nucleosome occupied (C2, C5, and C6), represented by aggregated V-plot NucleoATAC | Fig. 4c | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Nucleosome_phasing.ipynb) |
|  The six clusters included three types of V-plots where the central Etv2 sites were nucleosome free (C1, C3 and C4), and three types of V-plots where the central Etv2 sites were nucleosome occupied (C2, C5, and C6), represented by aggregated V-plot | Fig. 4d | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Nucleosome_phasing.ipynb) |
|  The six clusters included three types of V-plots where the central Etv2 sites were nucleosome free (C1, C3 and C4), and three types of V-plots where the central Etv2 sites were nucleosome occupied (C2, C5, and C6), represented by NucleoATAC | Fig. 4e | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Nucleosome_phasing.ipynb) |
|  The bar plot shows the proportion of each V-plot clusters in early, late and sustained Etv2 peaks, as well as in the background peaks | Fig. 4f | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Nucleosome_phasing.ipynb) |
|  The cluster labels of early, late and sustained ETV2 peak that were changed from MEF to D7 Flk1+ cells | Fig. 4g | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Nucleosome_phasing.ipynb) |
| The fragment size distribution of ATAC-seq of FLK1+ cells vs. the mixture population at 12 hours post-ETV2 induction during EB reprogramming and day 7 post-ETV2 induction during MEF reprogramming. In both conditions, the mono-nucleosomes and the di-nucleosomes were significantly increased in the FLK1+ cell populations | Fig. 4h | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/FLk1plus_mix.ipynb) |
| The aggregated V-plot whose centers are the ETV2 bound sites at FLK1+ cell populations at 12 hours post-ETV2 induction in EB, and at 7 days post-ETV2 induction | Fig. 4i | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/FLk1plus_mix.ipynb) |

### Figure 5

|  | Figures | Colab link | 
| --- | --- | --- | 
|  The UMAP plot shows the scRNA-seq of 8,838 cells from day 7 post-induction in control MEFs, 1,502 FLK1+ cells from day 7 post-ETV2 induction in MEFs, 8,248 cells from day 7 post-ETV2 induction in Brg1 KD MEFs, and 8,034 cells at day 7 in Brg1 KD MEFs  | Fig. 5c | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1KD_scRNA_seq_D7.ipynb) |
|  The UMAP plot shows the expression level of Etv2  | Fig. 5d | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1KD_scRNA_seq_D7.ipynb) |
|  The UMAP plot shows the expression level of Brg1  | Fig. 5e | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1KD_scRNA_seq_D7.ipynb) |
|  The UMAP plot shows the expression level of Flk1  | Fig. 5f | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1KD_scRNA_seq_D7.ipynb) |

### Figure 6

|  | Figures | Colab link | 
| --- | --- | --- | 
|  The heatmap shows the piled up ATAC-seq signal surrounding the summit of 12,170 sustained ETV2 ChIP-seq peaks that were present at day 1 and day 7 post-induction of ETV2 in control MEFs (the sustained Etv2 peaks) | Fig. 6a | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1_KD_sustained_Etv2_peaks.ipynb) |
|  The heatmap shows the transcription factors where motif associated chromatin accessibility were significantly changed at day 7 post-ETV2 induction in MEFs (unsorted MEFs or FLK1+ cells), or the Brg1 KD MEFs | Fig. 6b | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1_KD_sustained_Etv2_peaks.ipynb) |

### Figure 7

|  | Figures | Colab link | 
| --- | --- | --- | 
|  The heatmap shows chromatin accessibility for transcription factors for control vs. Brg1 knockout during ES/EB differentiation | Fig. 7f | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/chromVAR_analysis_of_Brg1_floxed_D4_and_Brg1_KO.ipynb) |
|  The heatmap shows the ChIP-seq of sustained ETV2 binding sites in EB differentiation, and ATAC-seq of Control and Brg1 KO EBs at day 4 | Fig. 7g | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Enriched_heatmap_of_Etv2_chip_seq_data_and_Brg1_floxed_Brg1_KO.ipynb) |
|  ATAC-seq V-plots of the genomic regions (640 bp) that are centered at ELK3 motifs | Fig. 7h | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Elk3_vplots.ipynb) |

## Extended Figures

### Extended Figure 1

|  | Figures | Colab link | 
| --- | --- | --- | 
|  The biological processes that are significantly associated with the up-regulated genes in cluster 7 (FLK1+ cells at day 7 of reprogramming) compared with cluster 1 (undifferentiated MEFs) | Extended Data Fig. 1h | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/scRNA_seq.ipynb) |

### Extended Figure 2

|  | Figures | Colab link | 
| --- | --- | --- | 
|  The violin plots show the scaled expression levels of endothelial markers such as Etv2, Emcn, Lmo2, Flk1/Kdr, Cdh5 and Sox18 in MEFs, day 1, day 2, day 7 post-ETV2 induction, as well as the FLK1+ cells from day 7 | Extended Data Fig. 2a | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/scRNA_seq.ipynb) |
|  The violin plots show the scaled expression levels of endothelial markers in seven cell clusters | Extended Data Fig. 2b | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_ChIP_seq_analysis.ipynb) |
|  The biological processes that are significantly associated with the up-regulated in genes in cluster 1 (undifferentiated MEFs) compared with the rest of the cell populations | Extended Data Fig. 2c | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/scRNA_seq.ipynb) |
|  GSEA plot indicates significant upregulation of the inflammatory response in MEFs (cluster 1) | Extended Data Fig. 2d | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/MEF_subpopulation_gsea_for_immune_reponse.ipynb) |
|  Heatmap representing the gene expression levels scaled by Seurat for upregulated (red) and downregulated (blue) genes in cluster 1 and cluster 2 | Extended Data Fig. 2e | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/R7_Etv2_pathway_MEF_clusters_1_and_2.ipynb) |
|  The bar plots show top 10 significant pathways for cluster 1 and cluster 2 for MEFs | Extended Data Fig. 2f | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/R7_Etv2_pathway_MEF_clusters_1_and_2.ipynb) |
|  The bar plots show top 10 significant pathways for cluster 1 and cluster 2 for MEFs. | Extended Data Fig. 2g | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/R7_Etv2_pathway_MEF_clusters_1_and_2.ipynb) |

### Extended Figure 3

|  | Figures | Colab link | 
| --- | --- | --- | 
|  The UMAP shows expression profiles of Tlr3, Nfkb1 and Vav3, Cd38 and Abl1 (members of B cell receptor signaling pathway) in undifferentiated MEFs and post Etv2 induction day 1, day2, day 7 and Flk1+ cells at day 7 | Extended Data Fig. 3a | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Expression_of_genes_related_Immune_reponse_pathways_in_the_7_clusters.ipynb) |
| The bar plot shows immune response related pathways significantly upregulated in Flk1+ cells from day 7 post Etv2 induction compared to MEFs | Extended Data Fig. 3b | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Gene_set_enrichment_analysis_Immune_reponse_pathways_for_MEFs_vs_D1_and_D7_Flk1%2B.ipynb) |
|  The bar plot shows immune response related pathways significantly upregulated in ay 1 post Etv2 induction compared to MEFs | Extended Data Fig. 3c | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Gene_set_enrichment_analysis_Immune_reponse_pathways_for_MEFs_vs_D1_and_D7_Flk1%2B.ipynb) |

### Extended Figure 5

|  | Figures | Colab link | 
| --- | --- | --- | 
|  The Venn diagrams show the overlap of commonly up-regulated genes during EB differentiation and MEF reprogramming | Extended Data Fig. 5a | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/scRNA_seq.ipynb) |
|  The Venn diagrams show the overlap of commonly down-regulated genes during EB differentiation and MEF reprogramming | Extended Data Fig. 5b | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/scRNA_seq.ipynb) |
|  Top commonly up-regulated genes during EB differentiation and MEF reprogramming | Extended Data Fig. 5c | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/scRNA_seq.ipynb) |
|  Top commonly down-regulated genes during EB differentiation and MEF reprogramming | Extended Data Fig. 5d | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/scRNA_seq.ipynb) |
|  The pathways that are significantly associated with commonly up-regulated genes during ES/EB differentiation and MEF reprogramming | Extended Data Fig. 5e | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/scRNA_seq.ipynb) |
|  The pathways that are significantly associated with commonly down-regulated genes during ES/EB differentiation and MEF reprogramming | Extended Data Fig. 5f | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/scRNA_seq.ipynb) |


### Extended Figure 6

|  | Figures | Colab link | 
| --- | --- | --- | 
|  The number of transcription factors whose motifs associated chromatin accessibility were significantly increased or decreased in the FLK1+ cell populations at 12 hours post-Etv2 induction compared with D2.5 EBs | Extended Data Fig. 6a | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/ATAC_analysis.ipynb) |
|  The number of transcription factors whose motifs associated chromatin accessibility were significantly increased or decreased in the FLK1+ cell population at day 7 post-ETV2 induction compared with undifferentiated MEFs | Extended Data Fig. 6b | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/ATAC_analysis.ipynb) |
|  The number of transcription factors whose motif associated chromatin accessibility that were commonly increased during EB and MEF reprogramming | Extended Data Fig. 6c | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/ATAC_analysis.ipynb) |
|  The number of transcription factors whose motif associated chromatin accessibility that were commonly decreased during EB and MEF reprogramming | Extended Data Fig. 6d | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/ATAC_analysis.ipynb) |
|  The transcription factors whose RNA-seq expression levels and motifs associated chromatin accessibility that were both up-regulated or down-regulated during EB reprogramming | Extended Data Fig. 6e | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/ATAC_analysis.ipynb) |
|  The transcription factors whose RNA-seq expression levels and motifs associated chromatin accessibility that were both up-regulated or down-regulated during MEF reprogramming | Extended Data Fig. 6f | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/ATAC_analysis.ipynb) |
|  The commonly up-regulated genes between EBs and MEFs | Extended Data Fig. 6g | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/TFBS_in_commonly_regulated_genes.ipynb) |
|  The commonly down-regulated genes between EBs and MEFs | Extended Data Fig. 6h | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/TFBS_in_commonly_regulated_genes.ipynb) |
|  The transcription factor motifs that are significantly enriched in 5k region surrounding the transcription start sites of the commonly up- and down-regulated genes in EBs and MEFs | Extended Data Fig. 6i | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/TFBS_in_commonly_regulated_genes.ipynb) |

### Extended Figure 7

|  | Figures | Colab link | 
| --- | --- | --- | 
|  The MNase-seq, BRG1, H3K27ac, H3, H3K9me3, H3K27me3, H3K9ac, H3K4me3, H4K7me1 and Hdac1 ChIP-seq signals surrounding the ETV2 bound sites at day1 post-ETV2 induction during MEF reprogramming, split into nucleosome and nucleosome free region (NFR) according to the MNase-seq signals | Extended Data Fig. 7a | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_binding_D1_MEF_extended.ipynb) |
|  The latent representation of ATAC-seq V-plots (-320bp to + 320bp) where the centers are nucleosome free or occupied by mono nucleosome | Extended Data Fig. 7b | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_ChIP_seq_analysis.ipynb) |
|  The aggregated ETV2 bound sites centric V-plot whose centers were occupied by mono nucleosomes or nucleosome free | Extended Data Fig. 7c | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_ChIP_seq_analysis.ipynb) |
|  The fragment size profile of ETV2 bound sites centric region (-320bp to +320bp) where the centers are nucleosome free or occupied by mono nucleosomes | Extended Data Fig. 7d | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_ChIP_seq_analysis.ipynb) |
|  Motif analysis of ETV2 bound sites during EB reprogramming and MEF reprogramming | Extended Data Fig. 7e | [R](diff_Etv2_motifs_between_NFR_and_nucleosome.Rmd) |
|  The overlap of ETV2 bound sites at day 1, day 2 and day 7 post-ETV2 induction during MEF reprogramming | Extended Data Fig. 7f | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_ChIP_seq_peaks.ipynb) |
|  The overlap of ETV2 bound sites at 3 hours and 12 hours post-ETV2 induction during EB reprogramming | Extended Data Fig. 7g | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_ChIP_seq_peaks.ipynb) |
|  The bar plot shows the percent of genes located near the late, early and sustained ETV2 bound sites related to blood vessel development | Extended Data Fig. 7h | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Pathway_Etv2_peaks.ipynb) |

### Extended Figure 9

|  | Figures | Colab link | 
| --- | --- | --- | 
|  The expression profile of Dek expression during MEF reprogramming | Extended Data Fig. 9a | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/Etv2_Project_Candidate_gene_Expression_Testing.ipynb) |
|  The expression profile of Znhit1 expression during MEF reprogramming | Extended Data Fig. 9b | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/Etv2_Project_Candidate_gene_Expression_Testing.ipynb) |
|  The expression profile of Chd8 expression during MEF reprogramming | Extended Data Fig. 9c | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/Etv2_Project_Candidate_gene_Expression_Testing.ipynb) |

### Extended Figure 10

|  | Figures | Colab link | 
| --- | --- | --- |
|  The heatmap shows the ATAC-seq signal surrounding the summit of 12,170 sustained ETV2 ChIP-seq peaks that were present at day 1 and day 7 post-induction of ETV2 in control MEFs | Extended Data Fig. 10a | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/) |
|  Heatmap shows ETV2 ChIP-seq signal surrounding 4,965 sustained ETV2 ChIP-seq peaks present in the D7 post ETV2 induction in WT MEFs and Brg1 KD MEFs | Extended Data Fig. 10b | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/Brg1_KD_ChIP_data_Enriched_Heatmap_at_Etv2_Binding_Sites.ipynb) |

## Supplementary Figure

### Supplementary Figure 1

|  | Figures | Colab link | 
| --- | --- | --- | 
|  The plots show enrichment scores for each Histone Acetyltransferases' | Supplementary Fig. 1a | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/R4_gsea_analysis_for_MEF_clusters.ipynb) |
|  The plots show enrichment scores for each Histone Deacetylases' | Supplementary Fig. 1b | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/R4_gsea_analysis_for_MEF_clusters.ipynb) |
|  The plots show enrichment scores for each Inflammatory Response | Supplementary Fig. 1c | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/R4_gsea_analysis_for_MEF_clusters.ipynb) |
|  The plots show enrichment scores for each NIK/NF-kappaB signaling | Supplementary Fig. 1d | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/R4_gsea_analysis_for_MEF_clusters.ipynb) |

### Supplementary Figure 2

|  | Figures | Colab link | 
| --- | --- | --- | 
|  The plots show enrichment scores for each Histone Acetyltransferases' | Supplementary Fig. 2a | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/R3ab_Copy_of_Expression_of_Inflammatory_related_genes_final.ipynb) <br> [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/R3_Gene_set_enrichment_analysis_MEF_vs_MEF_DOX_day1.ipynb)|
|  The plots show enrichment scores for each Histone Deacetylases' | Supplementary Fig. 2b | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/R3ab_Copy_of_Expression_of_Inflammatory_related_genes_final.ipynb) <br> [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/R3_Gene_set_enrichment_analysis_MEF_vs_MEF_DOX_day1.ipynb) |
|  The plots show enrichment scores for each Inflammatory Response | Supplementary Fig. 2c | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/R3ab_Copy_of_Expression_of_Inflammatory_related_genes_final.ipynb) <br> [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/R3_Gene_set_enrichment_analysis_MEF_vs_MEF_DOX_day1.ipynb) |
|  The plots show enrichment scores for each NIK/NF-kappaB signaling | Supplementary Fig. 2d | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/R3ab_Copy_of_Expression_of_Inflammatory_related_genes_final.ipynb) <br> [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/R3_Gene_set_enrichment_analysis_MEF_vs_MEF_DOX_day1.ipynb) |

### Supplementary Figure 3

|  | Figures | Colab link | 
| --- | --- | --- | 
|  The plots show enrichment scores for each Histone Acetyltransferases' | Supplementary Fig. 3a | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/R5_Expression_of_Inflammatory_pathway_related_genes_literature_data.ipynb) |
|  The plots show enrichment scores for each Histone Deacetylases' | Supplementary Fig. 3b | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/R5_Expression_of_Inflammatory_pathway_related_genes_literature_data.ipynb) |
|  The plots show enrichment scores for each Inflammatory Response | Supplementary Fig. 3c | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/R5_Expression_of_Inflammatory_pathway_related_genes_literature_data.ipynb) |
|  The plots show enrichment scores for each NIK/NF-kappaB signaling | Supplementary Fig. 3d | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/R5_Expression_of_Inflammatory_pathway_related_genes_literature_data.ipynb) |

### Supplementary Figure 4

|  | Figures | Colab link | 
| --- | --- | --- | 
|  Venn diagram representing the overlap of common genes upregulated in ES/EB and MEF reprogramming post Etv2 induction, also overlapping with the genes from gene ontology terms for Histone Acetyltransferases' | Supplementary Fig. 4a | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/Overlap_between_MEFs_and_EB_for_HDACs_HATs_Inflammatory_and_NFkB_genes.ipynb) |
|  Venn diagram representing the overlap of common genes upregulated in ES/EB and MEF reprogramming post Etv2 induction, also overlapping with the genes from gene ontology terms for Histone Deacetylases' | Supplementary Fig. 4b | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/Overlap_between_MEFs_and_EB_for_HDACs_HATs_Inflammatory_and_NFkB_genes.ipynb) |
|  Venn diagram representing the overlap of common genes upregulated in ES/EB and MEF reprogramming post Etv2 induction, also overlapping with the genes from gene ontology terms for Inflammatory Response | Supplementary Fig. 4c | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/Overlap_between_MEFs_and_EB_for_HDACs_HATs_Inflammatory_and_NFkB_genes.ipynb) |
|  Venn diagram representing the overlap of common genes upregulated in ES/EB and MEF reprogramming post Etv2 induction, also overlapping with the genes from gene ontology terms for NIK/NF-kappaB signaling | Supplementary Fig. 4d | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/Overlap_between_MEFs_and_EB_for_HDACs_HATs_Inflammatory_and_NFkB_genes.ipynb) |

### Supplementary Figure 5

|  | Figures | Colab link | 
| --- | --- | --- | 
|  The transcription factors whose motif associated chromatin accessibility and expressions were up-regulated in both EB and MEF on Etv2 induction | Supplementary Fig. 5a | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/ATAC_analysis.ipynb) |
|  The transcription factors whose motif associated chromatin accessibility and expressions were down-regulated in both EB and MEF on Etv2 induction | Supplementary Fig. 5b | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/ATAC_analysis.ipynb) |

### Supplementary Figure 6

|  | Figures | Colab link | 
| --- | --- | --- | 
|  The PCA analysis of undifferentiated MEFs and Brg1 KD MEFs | Supplementary Fig. 6a | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/ATAC_analysis.ipynb) |
|  The PCA analysis of different cell cycle programs | Supplementary Fig. 6b | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/ATAC_analysis.ipynb) |
|  The bar plot showing cells in the G1 phase are increased with Brg1 KD in MEFs | Supplementary Fig. 6c | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/ATAC_analysis.ipynb) |
|  The correction of cell cycle effects when combining the single cell RNA-seq data of undifferentiated MEFs and Brg1 KD MEFs | Supplementary Fig. 6d | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/ATAC_analysis.ipynb) |
|  The correction of cell cycle effects when combining the single cell RNA-seq data of undifferentiated MEFs and Brg1 KD MEFs | Supplementary Fig. 6e | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/ATAC_analysis.ipynb) |
|  The expression levels of Etv2 in cell-cycle effected adjusted single cell RNA-seq data from undifferentiated MEFs and Brg1 KD MEFs | Supplementary Fig. 6f | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/ATAC_analysis.ipynb) |
|  The expression levels of Etv2 in cell-cycle effected adjusted single cell RNA-seq data from undifferentiated MEFs and Brg1 KD MEFs | Supplementary Fig. 6g | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/ATAC_analysis.ipynb) |
|  The H3K27ac ChIP-seq signals of undifferentiated MEFs and Brg1 KD MEFs surrounding the early, late and sustained ETV2 bound sites in ETV2 induced MEF reprogramming | Supplementary Fig. 6h | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/H3K27ac_Chip_seq_Analysis.ipynb) |

## GEO datasets

The datasets below belong to a SuperSeries [GSE185684](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185684)
| Assay Type | GSE | 
| --- | --- |
| ATACseq| [GSE168636](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168636) | 
| ChIPseq | [GSE168521](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168521) | 
| RNA-seq |  [GSE185682](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185682) |  
| scRNA-seq | [GSE185683](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185683) |  
| NOMe-seq | [GSE185681](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185681) |  


### UCSC genome browser track

|  | Figures | UCSC Genome browser|
| --- | --- | --- |
| Lmo2 locus | Figure 2g | [Track](https://genome.ucsc.edu/s/gongx030/Etv2_pioneer_Lmo2) |
| All ChIP-seq data | | [Track](https://genome.ucsc.edu/s/gongx030%40umn.edu/Etv2_pioneer_ms) |
| Rhoj promoter | Figure 14c | [Track](https://genome.ucsc.edu/s/ndsouza/Etv2_Pioneer_ChIP_Rhoj) |
| Kdr promoter | Figure 14c | [Track](https://genome.ucsc.edu/s/ndsouza/Etv2_Pioneer_ChIP_Kdr) |


