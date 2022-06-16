#!/bin/sh
## From downloading RRBS data from GEO until Methylation base calling
##
##

cd data
mkdir rrbs
mkdir rrbs/raw_rrbs

## Downlaod RRBS data for the 25 samples
fastq-dump --skip-technical --outdir rrbs/raw_rrbs --gzip --split-files $(<RRBSSraAcclist.txt)

## QC check
for i in $(ls rrbs/raw_rrbs/*.fastq)
do
fastqc $i
done

## Trimming
cd rrbs
mkdir trimmed_rrbs

for i in $(ls raw_rrbs/*.fastq.gz)
do
fname="${i%%_*}"
trim_galore --rrbs -o trimmed_rrbs --paired ${fname}_1.fastq.gz ${fname}_2.fastq.gz
done

## Alignment
genome=~/UCSC/hg19/

cd rrbs
mkdir rrbs_bam

for i in $(ls trimmed_rrbs/*.fastq.gz)
fname="${i%%_*}"
do
bismark -o rrbs_bam $genome -1 ${fname}_1.fastq.gz -2 ${fname}_2.fastq.gz
done

## Methylation Base Calling
genome=~/UCSC/hg19/

for i in $(ls rrbs_bam)
do
bismark_methylation_extractor -p --comprehensive \
    --merge_non_CpG --bedGraph \
    --cytosine_report --genome_folder $genome \
    $i/${i}_1_bismark_bt2_pe.bam
done



