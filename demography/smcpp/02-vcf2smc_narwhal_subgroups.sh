#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evelien.degreef@umanitoba.ca
#SBATCH --mem=10GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%x-%j.out
#SBATCH --job-name=vcf2smc_run_narwhal_pops

source /home/edegreef/smcpp/bin/activate

cd /scratch/edegreef/narwhal_smcpp
mkdir masked_data_subgroup


list=scaf_min100kb

while read scaffold
do

# GF, RE, SB together (Canadian High Arctic, previously labeled WBB)
smc++ vcf2smc -d 94_RAHM_IQ_241_GF 94_RAHM_IQ_241_GF --mask smcpp_removed_sites_sorted.bed.gz narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz masked_data_subgroup/WBB.94_RAHM_IQ_241_GF.$scaffold.smc $scaffold WBB:94_RAHM_IQ_241_GF,94_RAHM_IQ_243_GF,94_RAHM_IQ_244_GF,ARGF_01_1068,ARGF_01_1082,ARGF_01_1088,ARGF_07_1121,ARGF_07_1125,ARGF_07_1127,ARGF_07_1131,ARRE_02_1094,ARRE_02_1106,ARRE_06_1144,ARRE_06_1146,ARRE_06_1157,ARRE_06_1158,ARRE_06_1159,ARRE_12_1295,ARSB_09_1050,ARSB_09_1053,B96_392_SB,B96_394_SB,B96_398_SB

# AB, CR, BI, PG, IG, PB group (Baffin Island, previously labeled EBB)
# mid/midblues
smc++ vcf2smc -d ARAB_04_03 ARAB_04_03 --mask smcpp_removed_sites_sorted.bed.gz narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz masked_data_subgroup/EBB.ARAB_04_03.$scaffold.smc $scaffold EBB:ARAB_04_03,ARAB_04_04,ARAB_08_1163,ARAB_08_1164,ARAB_10_1221,ARAB_10_1223,ARAB_99_1005,ARAB_99_1006,MM_AB_87_005,ARIG_06_1222,IMM_82_M1,ARPB_xx_1091,ARPB_xx_1095,ARPB_xx_1100,ARPG_19_1614,ARPI_04_1184,ARPI_08_0082,ARPI_17_1336,PI_92_009_MM,ARBI_05_1092,ARBI_05_1093,ARIQ_DFO_11_1071,94_RAHM_IQ_152,ARCR_07_1065,ARCR_07_1070

# Repulse Bay (Northern hudson bay/NHB)
smc++ vcf2smc -d ARRB_16_1379 ARRB_16_1379 --mask smcpp_removed_sites_sorted.bed.gz narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz masked_data_subgroup/RB.ARRB_16_1379.$scaffold.smc $scaffold RB:ARRB_16_1379,ARRB_99_1001,ARRB_99_1026,ARRB_xx_1421,ARRB_xx_1457,FMMM_RB_005_93,FMMM_RB_011_93,RBMM01,RBMM02

done < $list
