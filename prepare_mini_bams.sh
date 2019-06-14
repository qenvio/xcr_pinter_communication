#!/bin/bash

# set samtools version

alias samtools=/software/mb/el7.2/samtools-1.8/samtools

# function subset first Nth reads preserving header

function subset_bam {
	inbam=$1
	nreads=$2
	region=$3

	nhead=$(samtools view -H $inbam | wc -l)

	samtools view -h $inbam $region | head -n $(expr $nhead + $nreads) | samtools view -B
}



# Hi-C data

inbam=/users/project/4DGenome_no_backup/evidal/moritz/data/hic/bam/2019-01-21/ESH1_mask.bam
nreads=1000
region=X

subset_bam $inbam $nreads $region > hic_ESH1_mask_mini.bam

# ATACseq data

inbam=/users/project/4DGenome_no_backup/evidal/moritz/data/atacseq/bam/2018-04-27_ACCE61ANXX/ES1_mask.bam
nreads=1000
region=X

subset_bam $inbam $nreads $region > atacseq_ESH1_mask_mini.bam

# RNAseq data

inbam=/users/project/4DGenome_no_backup/evidal/moritz/data/rnaseq/bam/2018-08-10_ACCFWFANXX/ES1/ES1_all-clean.bam
nreads=1000
region=X

subset_bam $inbam $nreads $region > rnaseq_ESH1_mask_mini.bam
