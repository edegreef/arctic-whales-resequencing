#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evelien.degreef@umanitoba.ca
#SBATCH --mem=10GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%x-%j.out
#SBATCH --job-name=vcf2smc_run_bowredo

source /home/edegreef/smcpp/bin/activate

cd /scratch/edegreef/bowhead/smcpp_rightwhale
mkdir masked_data

list=scaf_min100kb_RW_autosomes

# running SMC with all samples as one population

while read scaffold
do
# Picking 2 different samples for distinguished lineages
smc++ vcf2smc -d 88_Pang 88_Pang --mask smcpp_removed_sites_sorted.bed.gz bowhead_RWmap_snps_rmfix.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz masked_data/ALL.88_Pang.$scaffold.smc $scaffold ALL:88_Pang,ARBMGH_2002_001,BMWG07_22,BMWG07_30,BM_01_2009,BM_CH_2000_01,BM_NSA_2008_02,BM_NSA_2009_02,BM_NSA_2009_03,BM_NSA_2010_01,BM_NSA_2011_01,BM_NSA_2011_03,BM_NSA_2012_02,BM_NSA_2012_03,BM_NSA_2014_01,BM_NSA_2020_01,BM_RB_2005_001,NSA_BM_98_01,WBF_2005_0298,BMDB_06_70

smc++ vcf2smc -d WBF_2005_0298 WBF_2005_0298 --mask smcpp_removed_sites_sorted.bed.gz bowhead_RWmap_snps_rmfix.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz masked_data/ALL.WBF_2005_0298.$scaffold.smc $scaffold ALL:88_Pang,ARBMGH_2002_001,BMWG07_22,BMWG07_30,BM_01_2009,BM_CH_2000_01,BM_NSA_2008_02,BM_NSA_2009_02,BM_NSA_2009_03,BM_NSA_2010_01,BM_NSA_2011_01,BM_NSA_2011_03,BM_NSA_2012_02,BM_NSA_2012_03,BM_NSA_2014_01,BM_NSA_2020_01,BM_RB_2005_001,NSA_BM_98_01,WBF_2005_0298,BMDB_06_70

done < $list
