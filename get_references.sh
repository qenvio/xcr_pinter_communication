#!/bin/bash

##  Settings

# prepare directory

basedir=~/moritz/data/reference/

mkdir -p $basedir/{fasta,vcf}


##  download data

# download VCFs

cd $basedir/vcf

wget -r --no-parent \
	 -nH --cut-dirs=3 \
	 -A '129S1_SvImJ.mgp.v5.*.dbSNP142*vcf.gz*' \
	 ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/

wget -r --no-parent \
	 -nH --cut-dirs=3 \
	 -A 'CAST_EiJ.mgp.v5.*.dbSNP142*vcf.gz*' \
	 ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/

wget -r --no-parent \
	 -nH --cut-dirs=3 \
	 -A 'mgp.v5.merged.*vcf.gz*' \
	 ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/


# download fastas

cd $basedir/fasta

wget -r --no-parent \
	 -nH --cut-dirs=3 \
	 -A '129S1_SvImJ.chromosomes.unplaced*' \
	 ftp://ftp-mouse.sanger.ac.uk/REL-1509-Assembly/

wget -r --no-parent \
	 -nH --cut-dirs=3 \
	 -A 'CAST_EiJ.chromosomes.unplaced*' \
	 ftp://ftp-mouse.sanger.ac.uk/REL-1509-Assembly/

wget -r --no-parent \
	 -nH --cut-dirs=3 \
	 -A 'GRCm38_68*' \
	 ftp://ftp-mouse.sanger.ac.uk/ref/

# merge both genomes in one fasta and index it

zcat 129S1_SvImJ.chromosomes.unplaced.gt2k.fa.gz |
	awk '{if(/^>/) $1 = $1"-129S1_SvImJ"; print}' > \
		both_strains.chromosomes.unplaced.gt2k.fa

zcat CAST_EiJ.chromosomes.unplaced.gt2k.fa.gz |
	awk '{if(/^>/) $1 = $1"-CAST_EiJ"; print}' >> \
		both_strains.chromosomes.unplaced.gt2k.fa

samtools faidx both_strains.chromosomes.unplaced.gt2k.fa

bwa index both_strains.chromosomes.unplaced.gt2k.fa

##  subset VCFs

bcftools view \
		 -O z -x -a \
		 -s 129S1_SvImJ,CAST_EiJ \
		 mgp.v5.merged.snps_all.dbSNP142.vcf.gz > \
		 mgp.v5.only_129S1_SvImJ-CAST_EiJ.snps_all.dbSNP142.vcf.gz

bcftools view \
		 -O z -x -a \
		 -s 129S1_SvImJ,CAST_EiJ \
		 mgp.v5.merged.indels.dbSNP142.normed.vcf.gz > \
		 mgp.v5.only_129S1_SvImJ-CAST_EiJ.indels.dbSNP142.normed.vcf.gz

# count SNPs

zgrep -v "^#" mgp.v5.only_129S1_SvImJ-CAST_EiJ.snps_all.dbSNP142.vcf.gz | \
	cut -f 10,11 | \
	awk '{gsub(":.*$", "", $1); gsub(":.*$", "", $2); print}' | \
	sort | \
	uniq -c | \
	sort -n | \
	sed -e 's,^ *\b,,g' -e 's, ,\t,g' > \
		mgp.v5.only_129S1_SvImJ-CAST_EiJ.snps_all.dbSNP142.counts.txt

sed -i 1i'n\t129S1_SvImJ\tCAST_EiJ' \
	mgp.v5.only_129S1_SvImJ-CAST_EiJ.snps_all.dbSNP142.counts.txt


## modify reference genome

# mask it with SNPs present in any of the two strains

bedtools maskfasta \
		 -fi ~/moritz/data/reference/fasta/GRCm38_68.fa \
		 -fo ~/moritz/data/reference/fasta/GRCm38_68_129S1_SvImJ-CAST_EiJ_snps_masked.fa \
		 -bed <(zgrep -v "^#" ~/moritz/data/reference/vcf/mgp.v5.only_129S1_SvImJ-CAST_EiJ.snps_all.dbSNP142.vcf.gz | awk -v OFS="\t" '{print $1, $2, $2}')

# index it

samtools faidx ~/moritz/data/reference/fasta/GRCm38_68_129S1_SvImJ-CAST_EiJ_snps_masked.fa

bwa index ~/moritz/data/reference/fasta/GRCm38_68_129S1_SvImJ-CAST_EiJ_snps_masked.fa


# unzip and index genomes

gunzip 129S1_SvImJ.chromosomes.unplaced.gt2k.fa.gz

samtools faidx 129S1_SvImJ.chromosomes.unplaced.gt2k.fa

gunzip CAST_EiJ.chromosomes.unplaced.gt2k.fa.gz

samtools faidx CAST_EiJ.chromosomes.unplaced.gt2k.fa.gz

# test


bcftools view -O z -r 19:10000000-11000000 \
		 mgp.v5.merged.snps_all.dbSNP142.vcf.gz > \
		 hola.vcf.gz

zgrep -cv "^#" hola.vcf.gz

bcftools view -s 129S1_SvImJ,CAST_EiJ hola.vcf.gz | zgrep -cv "^#"

bcftools view --trim-alt-alleles -s 129S1_SvImJ,CAST_EiJ hola.vcf.gz | zgrep -cv "^#"


##  repeated elements

mkdir -p $basedir/tracks

mysql -h genome-mysql.cse.ucsc.edu \
	  -u genome -D mm10 -N -A --column-names \
	  -e "select * from rmsk" |\
	gzip -c >\
		 $basedir/tracks/ucsc_mm10_rmsk.tsv.gz

zcat $basedir/tracks/ucsc_mm10_rmsk.tsv.gz |\
	awk -F "\t" -v OFS="\t" \
		'NR>1 && !($6 ~ /_/){gsub("chr", "", $6); print $6, $7, $8, $12, $13}' >\
		$basedir/tracks/ucsc_mm10_rmsk.bed

bgzip -c $basedir/tracks/ucsc_mm10_rmsk.bed >\
	  $basedir/tracks/ucsc_mm10_rmsk.bed.gz

tabix -p bed $basedir/tracks/ucsc_mm10_rmsk.bed.gz

## transcriptomes

masked_genome=~/moritz/data/reference/fasta/GRCm38_68_129S1_SvImJ-CAST_EiJ_snps_masked.fa

# get annotation from gencode

annot_dir=~/moritz/data/reference/annotation/
annot_file=$annot_dir/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf

mkdir -p $annot_dir
cd $annot_dir

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M19/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf.gz
gunzip gencode.vM19.chr_patch_hapl_scaff.annotation.gtf.gz


sed -i 's,^chr,,1' $annot_file

# index transcriptome with STAR

star_index_dir=~/moritz/data/reference/start/overhang_100
mkdir -p $star_index_dir


qsub -N star_index_masked_mm10\
	 -l virtual_free=32G\
	 -pe smp 8\
	 -j y\
	 -o $star_index_dir\
	 -b y\
	 star --runThreadN 8\
	 --runMode genomeGenerate\
	 --genomeDir $star_index_dir\
	 --genomeFastaFiles $masked_genome\
	 --sjdbGTFfile $annot_file\
	 --sjdbOverhang 100

ls $star_index_dir




