#!/bin/bash
## From downloading RNA-seq data from GEO until FPKM quantification
##
##

cd data
mkdir rna
mkdir rna/raw_rna

## Downlaod RNA-seq data for the 25 samples
fastq-dump --skip-technical --outdir rna/raw_rna --gzip --split-files $(<RNASraAcclist.txt)

## QC check
for i in $(ls rna/raw_rna/*.fastq)
do
fastqc $i
done

## Trimming 
mkdir rna/trimmed_rna

for i in $(ls raw_rna/*.fastq.gz)
do
fname="${i%%_*}"
java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar \
    PE ${fname}_1.fastq.gz ${fname}_2.fastq.gz \
    rna/trimmed_rna/${fname}_P1.fastq.gz untrimmed/${fname}_U1.fastq.gz \
    rna/trimmed_rna/${fname}_P2.fastq.gz untrimmed/${fname}_U2.fastq.gz \
    LEADING:20 SLIDINGWINDOW:5:20 \
    ILLUMINACLIP:~/Trimmomatic-0.39/adapters/nebnext.fa:2:30:10
done

## Alignment
genome=~/UCSC/hg19/genome.fa
gtf=~/UCSC/hg19/genes.gtf

cd rna
mkdir rna_bam

for i in $(ls raw_rna/*.fastq.gz)
do
fname="${i%%_*}"
tophat -p 8 -G $gtf -o rna_bam/${fname}_topout $genome ${fname}_P1.fastq.gz ${fname}_P2.fastq.gz
done

## Cufflinks
### Gene Expression Quantification and normalization

mkdir cuffq_out

for i in $(ls rna_bam/)
do
fname="${i%%_*}"
cuffquant --no-update-check -p 8 -o cuffq_out/${fname}_cqout $gtf ${i}/accepted_hits.bam
done

cuffnorm --no-update-check -p 8 --use-sample-sheet -o norm_out $gtf sample_sheet.txt