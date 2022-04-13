#!/bin/bash

# Map paired.fastq's to reference. This step takes a long time, example: ~5 days to make one 19GB bam file.
# Make sure reference genome is already indexed.

cd /home/degreefe/whales/BOW_fastq_trimmed/part6

# Going to do 6 parallel jobs at a time on biology-01
ls *_R1_paired.fastq.gz | sed 's/_R1_paired.fastq.gz$//' | parallel --jobs 6 'bwa mem /home/degreefe/whales/ref_genomes/BOW_reference.fasta {}_R1_paired.fastq.gz {}_R2_paired.fastq.gz | samtools sort -T {}_tmp -o {}.sorted.bam'

# Index bam file
ls *.sorted.bam | parallel 'samtools index {}'