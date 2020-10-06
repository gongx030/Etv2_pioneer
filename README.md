# The Etv2 Pioneer project


### [2020-10-02] `H3K27ac` analysis of Brg1 KD MEF

* Process the `H3K27ac` ChIP-seq data of Brg1 KD MEFs in [Alver et al.](https://www.nature.com/articles/ncomms14648#Sec6).  Use `Bowtie2` to map the reads to `mm10` and use `MACS2` to call the ChIP-seq peaks.  We will need the [`Fold Enrichment`](https://github.com/macs3-project/MACS/wiki/Build-Signal-Track) file from MACS2.

* `H3K27ac` ChIP-seq data of wildtype MEF: https://s3.msi.umn.edu/gongx030/datasets/dataset=Etv2PioneerChIPseq_version=20191203a/MEF_NoDox_d0_H3K27ac_FE.bw.  This bigwig file is the `Fold Enrichment` file from MACS2. 

* Compare with our `H3K27ac` ChIP-seq data in wildtype MEF at all Etv2 ChIP-seq binding sites.  The Etv2 ChIP-seq peaks can be found at https://s3.msi.umn.edu/gongx030/datasets/dataset=Etv2PioneerChIPseq_version=20191203a/all_Etv2_peaks.rds. This file contains aggregated Etv2 ChIP-seq peaks from ES/EB differentiation and MEF reprogramming.  Among all the peaks, we will be interested in the peaks that are present in D1 and D7 post Etv2 induction in MEF reprogramming.  Here is the code to get those peaks:

```
gr_url <- 'https://s3.msi.umn.edu/gongx030/datasets/dataset=Etv2PioneerChIPseq_version=20191203a/all_Etv2_peaks.rds' 
gr <- readRDS(gzcon(url(gr_url)))
peaks <- gr[gr$group[, 'MEF_Dox_d1_Etv2'] | gr$group[, 'MEF_Dox_d7_Etv2']]
peaks
```

## Processed datasets

| | link | Format | Script | 
| --- | --- | --- | --- | 
| ATAC-seq counts | [link](https://s3.msi.umn.edu/gongx030/datasets/dataset=Etv2ATAC_version=20190228a/all_ATAC.rds) | SummarizedExperiment | [ATAC_seq_preprocess.Rmd](ATAC_seq_preprocess.Rmd) |
| scRNA-seq of Etv2 induced reprogramming | [link](https://s3.msi.umn.edu/gongx030/etv2_pioneer/data/processed_Etv2_scRNAseq.rds) | SummarizedExperiment | [scRNA_seq_preprocess.Rmd](scRNA_seq_preprocess.Rmd) |
| Bulk RNA-seq of EB differentiation induced by Etv2 Dox | [link](https://s3.msi.umn.edu/gongx030/datasets/dataset=Etv2RNA-seq_version=20190909a/se.rds) | SummarizedExperiment | |

## Notebooks

|  | Figures | Preview | Colab link | Time |
| --- | --- | --- | --- | --- |
| Process scRNA-seq of<br> Etv2 reprogramming | | [R](scRNA_seq_preprocess.Rmd) |  | | |
| Process ATAC-seq data | | [R](ATAC_seq_preprocess.Rmd) |  |   |  |
| scRNA-seq of Etv2 induced reprogramming | Figure 1c-1g <br> Supplementary Figure 2 <br> Supplementray Figure 4 | [R](scRNA_seq.ipynb) | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/scRNA_seq.ipynb) | 11 mins |
| D7 Brg1 KD ATAC-seq | Figure 4a-4c | [R](Brg1_KD_sustained_Etv2_peaks.ipynb) | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1_KD_sustained_Etv2_peaks.ipynb) |  |
| D7 scRNA-seq including Brg1 KD samples | Figure 4e-4h | [R](Brg1KD_scRNA_seq_D7.ipynb) | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1KD_scRNA_seq_D7.ipynb) | 1.1 hrs |
| D0 scRNA-seq including Brg1 KD samples | Supplementary Figure 12 | [R](Brg1KD_scRNA_seq_D0.ipynb) | [R](https://colab.research.google.com/github/gongx030/etv2_pioneer/blob/master/Brg1KD_scRNA_seq_D0.ipynb) |  |


[ChIP_seq_preprocess.Rmd](ChIP_seq_preprocess.Rmd) R script for preprocessing Etv2 ChIP-seq, H3K27ac ChIP-seq and Brg1 ChIP-seq in ES/EB and MEF.  The combined ChIP-seq dataset sheet can be found [here](https://docs.google.com/spreadsheets/d/1UWiduM3Pv-GsVGmfxFApnyVBI1THMR8n8wHg5st3b5c/edit?usp=sharing).  

[diff_Etv2_motifs_between_NFR_and_nucleosome.Rmd](diff_Etv2_motifs_between_NFR_and_nucleosome.Rmd) R script for performing motif analysis of finding relatively enriched motifs in NFR compared with nucleosome within Etv2 ChIP-seq peaks regions. 

[find_de_novo_motifs_Etv2_chipseq_peaks.Rmd](find_de_novo_motifs_Etv2_chipseq_peaks.Rmd) R script for finding de novo motifs within Etv2 ChIP-seq peaks in EB and MEF.  The de novo motifs were identified by Homer. 

