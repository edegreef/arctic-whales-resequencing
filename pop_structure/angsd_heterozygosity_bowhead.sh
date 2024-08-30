#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evelien.degreef@umanitoba.ca
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=angsd_bw3
#SBATCH --output=%x-%j.out


module load StdEnv/2023 angsd/0.940

while read sample
do

angsd -i /scratch/edegreef/bowhead_bams/batch3/$sample.rightwhale.deDupRG.pp.bam -anc /scratch/edegreef/ref_genomes/GCA_028564815.2_mEubGla1.1.hap2._XY_genomic.fna -setMinDepthInd 5 -minmapq 25 -minq 25 -uniqueonly 1 -remove_bads 1 -GL 1 -dosaf 1 -rf scaf_min100kb_RW_autosomes -out $sample

realSFS $sample.saf.idx -fold 1 > $sample.est.ml

done < batch3