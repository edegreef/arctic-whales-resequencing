#!/bin/bash

#SBATCH --time=8-00:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evelien.degreef@umanitoba.ca
#SBATCH --constraint=skylake
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --job-name=map_fastq_part2
#SBATCH --output=%x-%j.out

# Load modules
module load StdEnv/2020 bwa/0.7.17 samtools/1.12

# Go to directory with merged & trimmed fastqs
cd /scratch/edegreef/whales/BOW_fastq/merged/trimmed

# Map paired.fastq's to reference. 
ls *_R1_paired.fastq.gz | sed 's/_R1_paired.fastq.gz$//' | parallel --jobs $SLURM_NTASKS 'bwa mem /scratch/edegreef/whales/ref_genomes/BOW_reference.fasta {}_R1_paired.fastq.gz {}_R2_paired.fastq.gz | samtools sort -T {}_tmp -o {}.sorted.bam'

# Index bam file
ls *.sorted.bam | parallel 'samtools index {}'
