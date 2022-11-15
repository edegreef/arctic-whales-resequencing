#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --cpus-per-task=12
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=downsample_batch1
#SBATCH --output=%x-%j.out

# This is part 1 of 2 (as in, doing this in 2 batches to speed it up)
# Downsampling some of the bam files to a lower coverage, aiming for 11x for the higher coverage samples (for bowhead)
# Proportion to 'keep' value based on target coverage divided by original coverage. for example if coverage=20x and want 10x, then keep=0.5

bam_ID1=88_Pang_S51
bam_ID2=AR_BM_SH_2003_01_S41
bam_ID3=ARBMGH_2002_001_S33
bam_ID4=BM_01_2009_S16
bam_ID5=BM_CH_2000_01_S32
bam_ID6=BM_NSA_2008_02_S38
bam_ID7=BM_NSA_2009_02_S50
bam_ID8=BM_NSA_2009_03_S36
bam_ID9=BM_NSA_2010_01_S39
bam_ID10=BM_NSA_2011_01_S42

keep1=0.6111
keep2=0.9167
keep3=0.9167
keep4=0.6875
keep5=0.6471
keep6=0.7333
keep7=0.7333
keep8=0.6471
keep9=0.6875
keep10=0.6111

cd /scratch/edegreef/whales/dedupRG_bam/samtools_filter
module load nixpkgs/16.09 gatk/4.1.2.0

gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID1.deDupRG.pp.bam -O downsampled/$bam_ID1.deDupRG.pp.downsampled.bam -P $keep1
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID2.deDupRG.pp.bam -O downsampled/$bam_ID2.deDupRG.pp.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID3.deDupRG.pp.bam -O downsampled/$bam_ID3.deDupRG.pp.downsampled.bam -P $keep3
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID4.deDupRG.pp.bam -O downsampled/$bam_ID4.deDupRG.pp.downsampled.bam -P $keep4
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID5.deDupRG.pp.bam -O downsampled/$bam_ID5.deDupRG.pp.downsampled.bam -P $keep5
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID6.deDupRG.pp.bam -O downsampled/$bam_ID6.deDupRG.pp.downsampled.bam -P $keep6
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID7.deDupRG.pp.bam -O downsampled/$bam_ID7.deDupRG.pp.downsampled.bam -P $keep7
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID8.deDupRG.pp.bam -O downsampled/$bam_ID8.deDupRG.pp.downsampled.bam -P $keep8
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID9.deDupRG.pp.bam -O downsampled/$bam_ID9.deDupRG.pp.downsampled.bam -P $keep9
gatk --java-options "-Xmx50G" DownsampleSam -I $bam_ID10.deDupRG.pp.bam -O downsampled/$bam_ID10.deDupRG.pp.downsampled.bam -P $keep10

# Check new modal coverage
module load StdEnv/2020 samtools/1.11

./bam_coverage.sh downsampled/$bam_ID1.deDupRG.pp.downsampled.bam
./bam_coverage.sh downsampled/$bam_ID2.deDupRG.pp.downsampled.bam
./bam_coverage.sh downsampled/$bam_ID3.deDupRG.pp.downsampled.bam
./bam_coverage.sh downsampled/$bam_ID4.deDupRG.pp.downsampled.bam
./bam_coverage.sh downsampled/$bam_ID5.deDupRG.pp.downsampled.bam
./bam_coverage.sh downsampled/$bam_ID6.deDupRG.pp.downsampled.bam
./bam_coverage.sh downsampled/$bam_ID7.deDupRG.pp.downsampled.bam
./bam_coverage.sh downsampled/$bam_ID8.deDupRG.pp.downsampled.bam
./bam_coverage.sh downsampled/$bam_ID9.deDupRG.pp.downsampled.bam
./bam_coverage.sh downsampled/$bam_ID10.deDupRG.pp.downsampled.bam

# Index all new bams at the end
module load StdEnv/2020 samtools/1.11

cd downsampled
for i in *.bam
do
echo "Indexing: "$i        
samtools index $i
done