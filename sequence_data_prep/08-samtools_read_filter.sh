#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=samtools_read_filter
#SBATCH --output=%x-%j.out

cd /scratch/edegreef/whales/BOW_bam/dedupRG

module load nixpkgs/16.09 gcc/8.3.0 samtools/1.9

# the -F 256 removes non-primary alignments. the -f 2 keeps only read mapped in proper pair

# one way, but accidentally turns files into "...deDupRG.bam.pp.bam"
for i in *deDupRG.bam
do
samtools view -b -h -F 256 -f 2 $i > samtools_filter/$i.pp.bam
done

