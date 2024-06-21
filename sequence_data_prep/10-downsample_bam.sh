#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evelien.degreef@umanitoba.ca
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=downsample_batch1
#SBATCH --output=%x-%j.out

# This script is for one downsampling batch of 6 samples. Ran this script with a few other batches simultaneously for other samples.

# List of sample name (including suffix from sequence lane)

bam_ID1=88_Pang_S51
bam_ID2=BMDB_06_70_S49
bam_ID3=BMWG07_22_S45
bam_ID4=BM_01_2009_S16_L001
bam_ID5=BM_CH_2000_01_S32
bam_ID6=BM_NSA_2008_02_S38
bam_ID7=ARBMGH_2002_001


# Keep1 corresponds to bam_ID1, keep2 corresponds to bam_ID2, etc...
# proportion of reads to keep (target coverage=25x)
keep1=0.6111
keep2=0.6875
keep3=0.6471
keep4=0.8462
keep5=0.6875
keep6=0.7333
keep7=0.92

# Go to directory with .deDupRG.pp.bam files  (also make sure the bam_coverage.sh script is in directory)
cd /scratch/edegreef/bowhead/rightwhale_map/bams_rightwhale/deDupRG/samtools_filter


# Load gatk
module load nixpkgs/16.09 gatk/4.1.2.0

# Run DownsampleSam command
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID1.rightwhale.deDupRG.pp.bam -O downsampled/$bam_ID1.rightwhale.deDupRG.pp.downsampled.bam -P $keep1
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID2.rightwhale.deDupRG.pp.bam -O downsampled/$bam_ID2.rightwhale.deDupRG.pp.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID3.rightwhale.deDupRG.pp.bam -O downsampled/$bam_ID3.rightwhale.deDupRG.pp.downsampled.bam -P $keep3
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID4.rightwhale.deDupRG.pp.bam -O downsampled/$bam_ID4.rightwhale.deDupRG.pp.downsampled.bam -P $keep4
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID5.rightwhale.deDupRG.pp.bam -O downsampled/$bam_ID5.rightwhale.deDupRG.pp.downsampled.bam -P $keep5
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID6.rightwhale.deDupRG.pp.bam -O downsampled/$bam_ID6.rightwhale.deDupRG.pp.downsampled.bam -P $keep6
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID7.rightwhale.deDupRG.pp.bam -O downsampled/$bam_ID7.rightwhale.deDupRG.pp.downsampled.bam -P $keep7



### Check new modal coverage
# Load samtools
module load StdEnv/2020 samtools/1.11

./bam_coverage.sh downsampled/$bam_ID1.rightwhale.deDupRG.pp.downsampled.bam
./bam_coverage.sh downsampled/$bam_ID2.rightwhale.deDupRG.pp.downsampled.bam
./bam_coverage.sh downsampled/$bam_ID3.rightwhale.deDupRG.pp.downsampled.bam
./bam_coverage.sh downsampled/$bam_ID4.rightwhale.deDupRG.pp.downsampled.bam
./bam_coverage.sh downsampled/$bam_ID5.rightwhale.deDupRG.pp.downsampled.bam
./bam_coverage.sh downsampled/$bam_ID6.rightwhale.deDupRG.pp.downsampled.bam
./bam_coverage.sh downsampled/$bam_ID7.rightwhale.deDupRG.pp.downsampled.bam
