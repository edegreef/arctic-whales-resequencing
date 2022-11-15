#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu
#SBATCH --constraint=skylake
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --job-name=deduplicate
#SBATCH --output=%x-%j.out

# Load modules
module load nixpkgs/16.09 picard/2.20.6

# Go to directory with sorted bam files
cd /scratch/edegreef/whales/BOW_bam

# Remove duplicate reads
ls *.sorted.bam | sed 's/.sorted.bam$//' | parallel --jobs 4 'java -Xmx50000m -jar $EBROOTPICARD/picard.jar MarkDuplicates I={}.sorted.bam O={}.deDup.bam M={}_deDupMetrics.txt REMOVE_DUPLICATES=true USE_JDK_DEFLATER=true USE_JDK_INFLATER=true'
