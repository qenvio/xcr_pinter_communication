#!/bin/bash

##  some "mock" commands with the alignment step
## (useless without the variables pintint to the files)

# ATACseq

bwa mem -t $ncores \
	$fa $fq1 $fq2 | \
	samtools view -b - | \
	samtools fixmate -@ $ncores -m - $outfile.tmp

# Hi-C

bwa mem -5SP -t $ncores \
	$fa $fq1 $fq2 | \
	samtools view -b - | \
	samtools fixmate -@ $ncores -m - $outfile.tmp

# RNAseq

star \
	--genomeDir $start_index_dir/ \
	--genomeLoad NoSharedMemory \
	--runThreadN $ncores \
	--outFilterType "BySJout" \
	--outFilterMultimapNmax 20 \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 \
	--outFilterMismatchNoverLmax 0.04 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--readFilesIn $fq1 $fq2 \
	--outSAMtype BAM SortedByCoordinate \
	--outTmpDir $outdir/tmp \
	--outFileNamePrefix $outdir/$id. \
	--outWigType bedGraph \
	--readFilesCommand zcat


##  versions

samtools --version
samtools 1.8
Using htslib 1.8
Copyright (C) 2018 Genome Research Ltd.


start --version
STAR_2.5.2b

bwa
Program: bwa (alignment via Burrows-Wheeler transformation)
Version: 0.7.15-r1142-dirty
Contact: Heng Li <lh3@sanger.ac.uk>
