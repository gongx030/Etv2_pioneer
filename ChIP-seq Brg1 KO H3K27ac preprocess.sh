#!/bin/bash
echo "installing dependencies"
apt install sra-toolkit
apt install bowtie2
echo "Importing Data"
wget https://genome-idx.s3.amazonaws.com/bt/mm10.zip #Get mm10 reference genome for mapping
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/chromInfo.txt.gz #get mm10 chromosome sizes
unzip mm10
gunzip chromInfo.txt.gz
fastq-dump SRR2133567
fastq-dump SRR2133571
mv SRR2133567.fastq SRR2133567_ChIP-seq_Smarca4_H3K27Acl.fastq #renamed to SRR2133567_ChIP-seq_Smarca4_H3K27Acl.fastq
mv SRR2133571.fastq SRR2133571_ChIP-seq_WT_H3K27AcI.fastq #renamed to SRR2133571_ChIP-seq_WT_H3K27AcI.fastq
echo "Mapping Reads"
for file in ./*.fastq
bowtie2 -t -x mm10 -U SRR2133567_ChIP-seq_Smarca4_H3K27Acl.fastq -S SRR2133567_ChIP-seq_Smarca4_H3K27Acl.sam -p 6 #CHANGE THREADS AS NEEDED
bowtie2 -t -x mm10 -U SRR2133571_ChIP-seq_WT_H3K27AcI.fastq -S SRR2133571_ChIP-seq_WT_H3K27AcI.sam -p 6 #CHANGE THREADS AS NEEDED
echo "converting to bam"
samtools view -S -b SRR2133567_ChIP-seq_Smarca4_H3K27Acl.sam > SRR2133567_ChIP-seq_Smarca4_H3K27Acl.bam
samtools view -S -b SRR2133571_ChIP-seq_WT_H3K27AcI.sam > SRR2133571_ChIP-seq_WT_H3K27AcI.bam
echo "clean and sort reads with picard" #(must install picard seperately)
picard SortSam I=SRR2133567_ChIP-seq_Smarca4_H3K27Acl.bam O=SRR2133567_ChIP-seq_Smarca4_H3K27Acl_Sorted.bam SORT_ORDER=coordinate
picard SortSam I=SRR2133571_ChIP-seq_WT_H3K27AcI.bam O=SRR2133571_ChIP-seq_WT_H3K27AcI_Sorted.bam SORT_ORDER=coordinate
picard FixMateInformation I=SRR2133567_ChIP-seq_Smarca4_H3K27Acl_Sorted.bam O=SRR2133567_ChIP-seq_Smarca4_H3K27Acl_Fixed_Mate.bam
picard FixMateInformation I=SRR2133571_ChIP-seq_WT_H3K27AcI_Sorted.bam O=SRR2133571_ChIP-seq_WT_H3K27AcI_Fixed_Mate.bam
picard MarkDuplicates I=RR2133567_ChIP-seq_Smarca4_H3K27Acl_Fixed_Mate.bam O=SRR2133567_ChIP-seq_Smarca4_H3K27Acl_markdup.bam M=SRR2133567_ChIP-seq_Smarca4_H3K27Acl_dupmetrics.txt REMOVE_DUPLICATES=false
picard MarkDuplicates I=SRR2133571_ChIP-seq_WT_H3K27AcI_Fixed_Mate.bam O=SRR2133571_ChIP-seq_WT_H3K27AcI_markdup.bam M=SRR2133571_ChIP-seq_WT_H3K27AcI_dupmetrics.txt REMOVE_DUPLICATES=false
picard BuildBamIndex I=SRR2133567_ChIP-seq_Smarca4_H3K27Acl_markdup.bam
picard BuildBamIndex I=SRR2133571_ChIP-seq_WT_H3K27AcI_markdup.bam
echo "Call Peaks"
macs2 callpeak -t RR2133567_ChIP-seq_Smarca4_H3K27Acl_markdup.bam -c SRR2133571_ChIP-seq_WT_H3K27AcI_markdup.bam -B --nomodel --broad --extsize 200 --shift 0 --keep-dup all -g mm -n ChIP-seq_Smarca4_H3K27Acl_FoldEnrichment
macs2 bdgcmp -t ChIP-seq_Smarca4_H3K27Acl_FoldEnrichment_treat_pileup.bdg -c ChIP-seq_Smarca4_H3K27Acl_FoldEnrichment_control_lambda.bdg -o ChIP-seq_Smarca4_H3K27Acl_FoldEnrichment_FE.bdg -m FE #Fold Enrichment
bdg2bw ChIP-seq_Smarca4_H3K27Acl_FoldEnrichment_FE.bdg chromInfo.txt #convert to bigwig
