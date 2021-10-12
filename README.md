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
| Differential Expressing testing candidate genes for knockdown in MEFs | | [R](https://colab.research.google.com/github/gongx030/Etv2_pioneer/blob/master/Etv2_Project_Candidate_gene_Expression_Testing.ipynb) |
| scRNA-seq of Etv2 induced reprogramming | Figure 1c-1g <br> Extended Data Fig.  2 <br> Extended Data Fig. 4 | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/scRNA_seq.ipynb) | 
| ATAC-seq and combined ATAC-seq/RNA-seq of <br> Etv2 induced MEF reprogramming and ES/EB differentiation | Figure 1h <br> Figure 1i <br> Extended Data Fig.  5 | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/ATAC_analysis.ipynb) | 
| Overlap of Etv2 ChIP-seq peaks in MEF and EB | Extended Data Fig.  8a and 8b | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_ChIP_seq_peaks.ipynb) | 
| Pathway analysis of early, late and sustained Etv2 peaks | Extended Data Fig.  8c | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Pathway_Etv2_peaks.ipynb) |
| Etv2 motifs in early, late and sustained Etv2 peaks | Figure 3a and 3b | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_motifs_in_early_Etv2_peaks.ipynb) |
| Early, late and sustained Etv2 peaks in **MEF** | Figure 3e | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/early_Etv2_peaks_in_MEF.ipynb) |
| D7 Brg1 KD ATAC-seq | Figure 4c | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1_KD_sustained_Etv2_peaks.ipynb) |
| chromVAR analysis of Brg1 KD ATAC-seq at D7| Figure 5a and 5b | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/chromVAR_Brg1_KD_ATAC_D7.ipynb) |
| D7 scRNA-seq including Brg1 KD samples | Figure 4c-4f | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1KD_scRNA_seq_D7.ipynb) | 
| D0 scRNA-seq including Brg1 KD samples | Extended Data Fig.  13a-13g | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1KD_scRNA_seq_D0.ipynb) |
| D0 ChIP-seq H3K27ac in Brg1 KD samples | Extended Data Fig.  13h | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/H3K27ac_Chip_seq_Analysis.ipynb) |
| Lmo2 binding sites | Figure 2g | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Lmo2_track.ipynb) |
| Enriched TFBS in commonly <br> up- and down-regulated genes <br> in EB and MEFs on Etv2 induction | Extended Data Fig.  5i, 5j, 5k | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/FBS_in_commonly_regulated_genes.ipynb) |
| NFR vs NOR at EB D2.5 | Figure 2d and 2b <br> Extended Data Fig.  6b-6d | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_ChIP_seq_analysis.ipynb) |  |  |
| NFR vs NOR at Etv2 binding sites at D1 MEF reprogramming | Figure 2c | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_binding_D1_MEF.ipynb) |
| NFR vs NOR at Etv2 binding sites at D1 MEF reprogramming <br> (including additional epigenetic marks) | Extended Data Fig.  6a | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Etv2_binding_D1_MEF_extended.ipynb) |
| D7 Brg1 KD ATAC-seq (NOR vs NFR) | Figure 4a | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1_KD_sustained_Etv2_peaks_NOR_NFR.ipynb) |
| scRNA-seq of Etv2 induced reprogramming at MEF and D1 | | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/scRNA_seq_D1.ipynb) | |
| We look at the enrichment of signals +/-1kb of Etv2 binding sites using Etv2 ChIP-seq 3h and 12h and Brg1 Control and KO ATAC-seq at Day 4. | Extended Data Fig.  15h | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Enriched_heatmap_of_Etv2_chip_seq_data_and_Brg1_floxed_Brg1_KO.ipynb) |
| chromVAR analysis for the change in chromatin accessibililty of Transcription factors at Brg1 Control and KO ATAC-seq at Day 4. | Extended Data Fig.  15g | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/chromVAR_analysis_of_Brg1_floxed_D4_and_Brg1_KO.ipynb) |
| NOMeseq analysis D0 and D1, MNase and Etv2 D1 enriched heatmap showinga few regions of the NOR cluster becoming NFR | Figure 2e| [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/2e_NOME_Enriched_Heatmap_with_MNase_and_Etv2.ipynb) |
| Enriched heatmap showing enrihment at early, late and sustained peaks for MEF Etv2 ChIP and Brg1 KD ChIP data | Extended Data Fig.  14b | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1_KD_ChIP_data_Enriched_Heatmap_at_Etv2_Binding_Sites.ipynb) |
| V-plots of Elk3 motif centric regions in D4 WT EB and D4 Brg1 KO EB | Figure 5f | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Elk3_vplots.ipynb) |
| Motif analysis of Etv2 ChIP-seq peaks | Figure 2f <br> Extended Data Fig. 7a | [R](find_de_novo_motifs_Etv2_chipseq_peaks.Rmd) [R](diff_Etv2_motifs_between_NFR_and_nucleosome.Rmd) |
| Dynamics of Etv2 peak centric V-plots during reprogramming | Extended Data Fig. 10d-10j | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Nucleosome_phasing.ipynb) |
| Fragment size distribution between Flk1+ cells and unsorted mixture cells | Extended Data Fig. 10a-10c | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/FLk1plus_mix.ipynb) |

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

