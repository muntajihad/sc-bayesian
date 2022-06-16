#!/bin/bash
## Estimationg CNV from RRBS data using HMMcopy
## Part1: Generating wig files.
##
##
mkdir cnv
mkdir cnv/wigs
## Generate wig file for genomic gc content and mappability.
genome=~/UCSC/hg19/genome.fa
map=~/UCSC/hg19/wgEncodeCrgMapabilityAlign100mer.bigWig

mapCounter -w 100000 $map > cnv/wigs/map.wig
gcCounter -w 100000 $genome > cnv/wigs/gc.wig

## In order to generate wig files, BAM files should be sorted and indexed first.
cd rrbs/rrbs_bam
mkdir sorted

for f in *.bam
do
fname="${f%%.*}"
samtools sort -@ 7 $f sorted/${fname}.sorted.bam
samtools index sorted/${fname}.sorted.bam
done

## Generate wig files for each sample from sorted BAM files
cd sorted

for f in *.bam
do
filename="${f%%.*}"
readCounter -w 100000 $f > ~/cnv/wigs/${filename}.wig
done

# Remove wingows from mito DNA 
cd cnv/wigs

for i in $(ls )
do
tail -n +5 $i | sponge $i
done
