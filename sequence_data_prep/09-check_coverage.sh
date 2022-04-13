#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=bam_coverage
#SBATCH --output=%x-%j.out

# Load modules
module load StdEnv/2020 samtools/1.12

# Go to directory with final bam files
cd /scratch/edegreef/whales/BOW_bam/dedupRG

# Check modal coverage for all deDupRG.bam files in folder, using code from Phil Grayson for the bam_coverage.sh script (https://github.com/phil-grayson/SexFindR/blob/main/Step_1/samtools_Modaldepth.sh)
for i in *deDupRG.bam
do
echo "Running bam_coverage: "$i        
./bam_coverage.sh $i
done

# Can also do one file at a time (example)
#./bam_coverage.sh Ha14-13-716-505.deDupRG.bam
