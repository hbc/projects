#!/bin/bash
for file in *522*.bam
do
    stem=`basename $file .bam`
    echo "Intersecting $file with rRNA genes."
    bedtools intersect -abam $file -wa -b ref-transcripts-rRNA.gtf > region/${stem}_rRNA.bam
    samtools view -c region/${stem}_rRNA.bam > region/${stem}_rRNA.count
    echo "Intersecting $file with PD genes."
    bedtools intersect -abam $file -wa -b ref-transcripts-PD.gtf > region/${stem}_PD.bam
    samtools view -c region/${stem}_PD.bam > region/${stem}_PD.count
    echo "Intersecting $file with transcriptome."
    bedtools intersect -abam $file -wa -b ref-transcripts.gtf > region/${stem}_transcriptome.bam
    samtools view -c region/${stem}_transcriptome.bam > region/${stem}_transcriptome.count
    echo "Counting reads mapping to transcriptome in $file."
    samtools view -h region/${stem}_transcriptome.bam > region/${stem}_transcriptome.sam
    htseq-count --stranded=no region/${stem}_transcriptome.sam ref-transcripts.gtf > region/counts/${stem}_transcriptome.counts
 done
