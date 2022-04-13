#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu
#SBATCH --constraint=skylake
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --job-name=addRG
#SBATCH --output=%x-%j.out

# Load modules
module load nixpkgs/16.09 gcc/8.3.0 picard/2.20.6 samtools/1.9

# Create reference genome dictionary (.dict file)
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=/scratch/edegreef/ref_genomes/BOW_reference.fasta O=/scratch/edegreef/ref_genomes/BOW_reference.dict 

# Add read group information to each deDup bam
ls /scratch/edegreef/whales/BOW_bam/*deDup.bam | sed 's/.deDup.bam$//' | parallel --jobs 8 'java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I={}.deDup.bam O={}.deDupRG.bam RGID={} RGPL=illumina RGSM={} RGLB=lib1 RGPU=unit1'

# Index new bam files
ls /scratch/edegreef/whales/BOW_bam/*deDupRG.bam | sed 's/.deDupRG.bam$//' | parallel --jobs 8 'samtools index {}.deDupRG.bam'