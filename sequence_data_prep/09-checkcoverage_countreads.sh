#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu
#SBATCH --mem=35G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=bam_coverage
#SBATCH --output=%x-%j.out

# Load modules
module load StdEnv/2020 samtools/1.12

# Go to directory with final bam files
cd /scratch/edegreef/narwhal/bams_sorted

# Check modal coverage for all deDupRG.bam files in folder, using code from Phil Grayson for the bam_coverage.sh script (https://github.com/phil-grayson/SexFindR/blob/main/Step_1/samtools_Modaldepth.sh)
for i in *sorted.bam
do
echo "Running bam_coverage: "$i        
./bam_coverage.sh $i

# also count reads
samtools view -c $i > $i.ctotalreads
#samtools view -F 0x04 -c $i > $i.cmappedreads

samtools depth $i  |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > $i.samtools.depth

done

# Can also do one file at a time (example)
#./bam_coverage.sh Ha14-13-716-505.deDupRG.bam
