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

## Datasets

The processed sequencing data can be found [here](https://docs.google.com/spreadsheets/d/18E90G_y3H5jT5wbNObcNY0Jdhqo0YWRMrPHzvoaiwD8/edit?usp=sharing). 

## Notebooks

[Brg1KD_scRNA-seq.ipynb](https://colab.research.google.com/drive/18I1xy3uPzhEfTzZZbM-bPJTITMnPoHHB?usp=sharing) Analysis of Brg1KD single cell RNA-seq. 

[Brg1KD_analysis.ipynb](Brg1KD_analysis.ipynb) Notebook for Brg1 KD ATAC-seq analysis. 

[ATAC_seq_preprocess.Rmd](ATAC_seq_preprocess.Rmd) R script for processing ATAC-seq data of ES/EB differentiation and MEF reprogramming induced by Etv2.  

[ChIP_seq_preprocess.Rmd](ChIP_seq_preprocess.Rmd) R script for preprocessing Etv2 ChIP-seq, H3K27ac ChIP-seq and Brg1 ChIP-seq in ES/EB and MEF.  The combined ChIP-seq dataset sheet can be found [here](https://docs.google.com/spreadsheets/d/1UWiduM3Pv-GsVGmfxFApnyVBI1THMR8n8wHg5st3b5c/edit?usp=sharing).  

[scRNA_seq_analysis.Rmd](scRNA_seq_analysis.Rmd) R script for performing the scRNA-seq of Etv2 reprogramming in MEF. 

[diff_Etv2_motifs_between_NFR_and_nucleosome.Rmd](diff_Etv2_motifs_between_NFR_and_nucleosome.Rmd) R script for performing motif analysis of finding relatively enriched motifs in NFR compared with nucleosome within Etv2 ChIP-seq peaks regions. 

[find_de_novo_motifs_Etv2_chipseq_peaks.Rmd](find_de_novo_motifs_Etv2_chipseq_peaks.Rmd) R script for finding de novo motifs within Etv2 ChIP-seq peaks in EB and MEF.  The de novo motifs were identified by Homer. 
