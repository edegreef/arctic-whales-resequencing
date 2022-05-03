#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@myucdavis.edu
#SBATCH --mem=10GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%x-%j.out
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=vcf2smc_run1

source /home/edegreef/smcpp/bin/activate

cd /scratch/edegreef/whales/BOW_smcpp
mkdir masked_data_run_1group

list=scaf_min100kb

# running SMC with all samples as one population. For now just using first sample in list as distinguished lineage

while read scaffold
do
# Labrador as main
smc++ vcf2smc -d 88_Pang 88_Pang --mask smcpp_removed_sites_sorted.bed.gz BOW_SNPS.filter1.miss.biallel.min100kb.autosomes.vcf.gz masked_data_run_1group/ALL.88_Pang.$scaffold.smc $scaffold ALL:88_Pang,99_01,AR_BM_SH_2003_01,ARBMGH_2002_001,BM_01_2009,BM_CH_2000_01,BM_NSA_2008_02,BM_NSA_2009_02,BM_NSA_2009_03,BM_NSA_2010_01,BM_NSA_2011_01,BM_NSA_2011_03,BM_NSA_2012_02,BM_NSA_2012_03,BM_NSA_2014_01,BM_NSA_2020_01,BM_RB_2005_001,BMDB_06_70,BMWG07_22,BMWG07_30,NSA_BM_98_01,RMD_BM_96_1,WBF_2005_0298

done < $list
