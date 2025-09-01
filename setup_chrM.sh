#!/usr/bin/env bash
set -e
mkdir -p data
cd data
aria2c -x 8 -s 8 -o chrM.fa.gz https://hgdownload.soe.ucsc.edu/goldenpath/hg38/chromosomes/chrM.fa.gz
gunzip -f chrM.fa.gz
samtools faidx chrM.fa
bwa index chrM.fa
wgsim -N 5000 -1 100 -2 100 chrM.fa sample_R1.fastq sample_R2.fastq
gzip -f sample_R1.fastq sample_R2.fastq
echo "Reference + reads ready"
