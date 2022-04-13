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
#SBATCH --job-name=platypus
#SBATCH --output=%x-%j.out

# Load modules
module load nixpkgs/16.09 gcc/7.3.0 platypus/0.8.1

# Go to directory with final bam files
cd /scratch/edegreef/whales/BOW_bam/dedupRG

# Run platypus to call variants from bams & reference genome 
# Make sure all bam files are indexed beforehand. Also couldn't add Platypus.py to path for some reason so I inputed the path here in script
python /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/gcc7.3/platypus/0.8.1/bin/Platypus.py callVariants \
--bamFiles=88_Pang_S51.deDupRG.bam,99_01_S15.deDupRG.bam,AR_BM_SH_2003_01_S41.deDupRG.bam,ARBMGH_2002_001_S33.deDupRG.bam,BM_01_2009_S16.deDupRG.bam,BM_CH_2000_01_S32.deDupRG.bam,BM_NSA_2008_02_S38.deDupRG.bam,BM_NSA_2009_02_S50.deDupRG.bam,BM_NSA_2009_03_S36.deDupRG.bam,BM_NSA_2010_01_S39.deDupRG.bam,BM_NSA_2011_01_S42.deDupRG.bam,BM_NSA_2011_03_S44.deDupRG.bam,BM_NSA_2012_02_S47.deDupRG.bam,BM_NSA_2012_03_S34.deDupRG.bam,BM_NSA_2014_01_S17.deDupRG.bam,BM_NSA_2020_01_S35.deDupRG.bam,BM_RB_2005_001_S43.deDupRG.bam,BMDB_06_70_S49.deDupRG.bam,BMWG07_22_S45.deDupRG.bam,BMWG07_30_S48.deDupRG.bam,NSA_BM_98_01_S40.deDupRG.bam,RMD_BM_96_1_S37.deDupRG.bam,WBF_2005_0298_S46.deDupRG.bam \
--refFile=/scratch/edegreef/whales/ref_genomes/BOW_reference.fasta \
--output=/scratch/edegreef/whales/BOW_snps/BOW_allvariantcalls.vcf \
--nCPU=48