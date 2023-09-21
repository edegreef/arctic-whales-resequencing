#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evelien.degreef@umanitoba.ca
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --nodes=1
#SBATCH --job-name=platypus
#SBATCH --output=%x-%j.out

# Load modules
module load nixpkgs/16.09 gcc/7.3.0 platypus/0.8.1

# set up for bowhead whale
# same set up for narwhal except for using "narwhal_bams.txt" and "NAR_reference.fasta" and output for narwhal name.

# Go to directory with final bam files (in hindsight can just make path to each bam file in the .txt list file instead of going to folder).
cd /scratch/edegreef/whales/dedupRG_bam/samtools_filter/downsampled

# Run platypus to call variants from bams & reference genome 
# Make sure all bam files are indexed beforehand. Also couldn't add Platypus.py to path for some reason so I inputed the path here in script
python /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/gcc7.3/platypus/0.8.1/bin/Platypus.py callVariants \
--bamFiles=bowhead_whale_bams.txt \
--refFile=/scratch/edegreef/whales/ref_genomes/rclone/BOW_reference.fasta \
--output=/scratch/edegreef/whales/bowhead_allvariantcalls_pp.vcf \
--nCPU=32 --minReads=4
