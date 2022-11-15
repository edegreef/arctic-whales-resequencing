#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --nodes=1
#SBATCH --job-name=platypus
#SBATCH --output=%x-%j.out

# Load modules
module load nixpkgs/16.09 gcc/7.3.0 platypus/0.8.1

# Go to directory with final bam files
cd /scratch/edegreef/whales/dedupRG_bam/samtools_filter/downsampled

# Run platypus to call variants from bams & reference genome 
# Make sure all bam files are indexed beforehand. Also couldn't add Platypus.py to path for some reason so I inputed the path here in script
python /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/gcc7.3/platypus/0.8.1/bin/Platypus.py callVariants \
--bamFiles=88_Pang_S51.deDupRG.pp.downsampled.bam,/scratch/edegreef/whales/dedupRG_bam/samtools_filter/99_01_S15.deDupRG.pp.bam,ARBMGH_2002_001_S33.deDupRG.pp.downsampled.bam,BM_01_2009_S16.deDupRG.pp.downsampled.bam,AR_BM_SH_2003_01_S41.deDupRG.pp.downsampled.bam,BM_CH_2000_01_S32.deDupRG.pp.downsampled.bam,BMDB_06_70_S49.deDupRG.pp.downsampled.bam,BM_NSA_2008_02_S38.deDupRG.pp.downsampled.bam,BM_NSA_2009_02_S50.deDupRG.pp.downsampled.bam,BM_NSA_2009_03_S36.deDupRG.pp.downsampled.bam,BM_NSA_2010_01_S39.deDupRG.pp.downsampled.bam,BM_NSA_2011_01_S42.deDupRG.pp.downsampled.bam,BM_NSA_2011_03_S44.deDupRG.pp.downsampled.bam,BM_NSA_2012_02_S47.deDupRG.pp.downsampled.bam,BM_NSA_2012_03_S34.deDupRG.pp.downsampled.bam,/scratch/edegreef/whales/dedupRG_bam/samtools_filter/BM_NSA_2014_01_S17.deDupRG.pp.bam,BM_NSA_2020_01_S35.deDupRG.pp.downsampled.bam,BM_RB_2005_001_S43.deDupRG.pp.downsampled.bam,BMWG07_22_S45.deDupRG.pp.downsampled.bam,BMWG07_30_S48.deDupRG.pp.downsampled.bam,/scratch/edegreef/whales/dedupRG_bam/samtools_filter/NSA_BM_98_01_S40.deDupRG.pp.bam,RMD_BM_96_1_S37.deDupRG.pp.downsampled.bam,WBF_2005_0298_S46.deDupRG.pp.downsampled.bam \
--refFile=/scratch/edegreef/whales/ref_genomes/rclone/BOW_reference.fasta \
--output=/scratch/edegreef/whales/bowhead_allvariantcalls_pp.vcf \
--nCPU=32 --minReads=4
