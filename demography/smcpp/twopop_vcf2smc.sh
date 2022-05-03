#!/bin/bash

#SBATCH --time=24:00:00 
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@myucdavis.edu
#SBATCH --mem=10GB 
#SBATCH --ntasks=1
#SBATCH --output=%x-%j.out
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=vcf2smc_run_part3

source /home/edegreef/smcpp/bin/activate

cd /scratch/edegreef/whales/NAR_smcpp

mkdir masked_data_run_ABGF

list=scaf_min100kb

# running SMC with pairwise pops

while read scaffold
do

# Pop1=AB, Pop2=GF
smc++ vcf2smc -d ARAB_04_03 ARAB_04_03 --mask smcpp_removed_sites_sorted.bed.gz NAR_SNPS.filter1.miss.biallel.min100kb.autosomes.vcf.gz masked_data_run_ABGF/ABGF.$scaffold.smc $scaffold AB:ARAB_04_03,ARAB_04_04,ARAB_08_1163,ARAB_08_1164,ARAB_10_1221,ARAB_10_1223,ARAB_99_1005,ARAB_99_1006,MM_AB_87_005 GF:94_RAHM_IQ_241_GF,94_RAHM_IQ_243_GF,94_RAHM_IQ_244_GF,ARGF_01_1068,ARGF_01_1082,ARGF_01_1088,ARGF_07_1121,ARGF_07_1125,ARGF_07_1127,ARGF_07_1131


done < $list