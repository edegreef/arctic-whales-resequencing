#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@myucdavis.edu
#SBATCH --output=%x-%j.out

# load modules
module load StdEnv/2020 fastqc/0.11.9

# go to directory with the raw fastqs
cd BOW_fastq

# run fastqc, using 6 threads here
fastqc -t 6 *fastq.gz

