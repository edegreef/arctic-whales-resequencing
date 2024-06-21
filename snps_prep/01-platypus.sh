#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evelien.degreef@umanitoba.ca
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --nodes=1
#SBATCH --job-name=platypus
#SBATCH --output=%x-%j.out

# Go to directory with final bam files
cd /scratch/edegreef/bowhead/rightwhale_map

# Load modules
module load nixpkgs/16.09 gcc/7.3.0 platypus/0.8.1

# Run platypus to call variants from bams & reference genome 
# Make sure all bam files are indexed beforehand. Also couldn't add Platypus.py to path for some reason so I inputed the path here in script
python /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/gcc7.3/platypus/0.8.1/bin/Platypus.py callVariants \
--bamFiles=bowhead_rightwhale_bams_list.txt \
--refFile=/scratch/edegreef/bowhead/ref_genome/GCA_028564815.2_mEubGla1.1.hap2._XY_genomic.fna \
--output=/scratch/edegreef/bowhead/rightwhale_map/bowhead_RWmap_allvariants.vcf \
--nCPU=32 --minReads=4

