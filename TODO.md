
### [2020-10-02] `H3K27ac` analysis of Brg1 KD MEF
================================================================================
* Process the `H3K27ac` ChIP-seq data of Brg1 KD MEFs in [Alver et al.](https://www.nature.com/articles/ncomms14648#Sec6).  Use `Bowtie2` to map the reads to `mm10` and use `MACS2` to call the ChIP-seq peaks.  We will need the [`Fold Enrichment`](https://github.com/macs3-project/MACS/wiki/Build-Signal-Track) file from MACS2.

* `H3K27ac` ChIP-seq data of wildtype MEF: https://s3.msi.umn.edu/gongx030/datasets/dataset=Etv2PioneerChIPseq_version=20191203a/MEF_NoDox_d0_H3K27ac_FE.bw.  This bigwig file is the `Fold Enrichment` file from MACS2. 

* Compare with our `H3K27ac` ChIP-seq data in wildtype MEF at all Etv2 ChIP-seq binding sites.  The Etv2 ChIP-seq peaks can be found at https://s3.msi.umn.edu/gongx030/datasets/dataset=Etv2PioneerChIPseq_version=20191203a/all_Etv2_peaks.rds. This file contains aggregated Etv2 ChIP-seq peaks from ES/EB differentiation and MEF reprogramming.  Among all the peaks, we will be interested in the peaks that are present in D1 and D7 post Etv2 induction in MEF reprogramming.  Here is the code to get those peaks:

```
gr_url <- 'https://s3.msi.umn.edu/gongx030/datasets/dataset=Etv2PioneerChIPseq_version=20191203a/all_Etv2_peaks.rds' 
gr <- readRDS(gzcon(url(gr_url)))
peaks <- gr[gr$group[, 'MEF_Dox_d1_Etv2'] | gr$group[, 'MEF_Dox_d7_Etv2']]
peaks
```
