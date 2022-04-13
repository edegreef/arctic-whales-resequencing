#!/bin/bash

# Using bedtools to prepare reference bed file and annotate which windows are marked sex-linked and control

# Input reference genome fai path and ref genome prefix
ref_genome_path=/home/degreefe/whales/ref_genomes/
ref_genome=BOW_reference

# Make bed file for reference genome using the fasta.fai file
bedtools makewindows -g $ref_genome_path/$ref_genome.fasta.fai -w 10000 > $ref_genome.10kb.bed

# Run annotate for X scaffolds using bed files for each combo M1,M2,F1,F2 for X-linked ones
bedtools annotate -i $ref_genome.10kb.bed -files X_M1F1.bed X_M1F2.bed X_M2F1.bed X_M2F2.bed control_M1M2.bed control_F1F2.bed > $ref_genome.10kb.Xlinked.bed

# Run annotate for Y scaffolds using bed files for each combo M1,M2,F1,F2 for Y-linked ones
bedtools annotate -i $ref_genome.10kb.bed -files Y_M1F1.bed Y_M1F2.bed Y_M2F1.bed Y_M2F2.bed control_M1M2.bed control_F1F2.bed > $ref_genome.10kb.Ylinked.bed