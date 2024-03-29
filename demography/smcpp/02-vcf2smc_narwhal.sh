#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evelien.degreef@umanitoba.ca
#SBATCH --mem=10GB
#SBATCH --ntasks=1
#SBATCH --output=%x-%j.out
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=vcf2smc_run_narredo

source /home/edegreef/smcpp/bin/activate

cd /scratch/edegreef/whales/NAR_smcpp
mkdir masked_data

list=scaf_min100kb

while read scaffold
do
# Running on 2 different distinguished lineages
smc++ vcf2smc -d 94_RAHM_IQ_152 94_RAHM_IQ_152 --mask smcpp_removed_sites_sorted.bed.gz narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz masked_data/ALL.94_RAHM_IQ_152.$scaffold.smc $scaffold ALL:94_RAHM_IQ_152,94_RAHM_IQ_241_GF,94_RAHM_IQ_243_GF,94_RAHM_IQ_244_GF,ARAB_04_03,ARAB_04_04,ARAB_08_1163,ARAB_08_1164,ARAB_10_1221,ARAB_10_1223,ARAB_99_1005,ARAB_99_1006,ARBI_05_1092,ARBI_05_1093,ARCR_07_1065,ARCR_07_1070,ARGF_01_1068,ARGF_01_1082,ARGF_01_1088,ARGF_07_1121,ARGF_07_1125,ARGF_07_1127,ARGF_07_1131,ARIG_06_1222,ARIQ_DFO_11_1071,ARPB_xx_1091,ARPB_xx_1095,ARPB_xx_1100,ARPG_19_1614,ARPI_04_1184,ARPI_08_0082,ARPI_17_1336,ARRB_16_1379,ARRB_99_1001,ARRB_99_1026,ARRB_xx_1421,ARRB_xx_1457,ARRE_02_1094,ARRE_02_1106,ARRE_06_1144,ARRE_06_1146,ARRE_06_1157,ARRE_06_1158,ARRE_06_1159,ARRE_12_1295,ARSB_09_1050,ARSB_09_1053,B96_392_SB,B96_394_SB,B96_398_SB,FMMM_RB_005_93,FMMM_RB_011_93,IMM_82_M1,MM_AB_87_005,PI_92_009_MM,RBMM01,RBMM02


smc++ vcf2smc -d ARGF_01_1068 ARGF_01_1068 --mask smcpp_removed_sites_sorted.bed.gz narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz masked_data/ALL.ARGF_01_1068.$scaffold.smc $scaffold ALL:94_RAHM_IQ_152,94_RAHM_IQ_241_GF,94_RAHM_IQ_243_GF,94_RAHM_IQ_244_GF,ARAB_04_03,ARAB_04_04,ARAB_08_1163,ARAB_08_1164,ARAB_10_1221,ARAB_10_1223,ARAB_99_1005,ARAB_99_1006,ARBI_05_1092,ARBI_05_1093,ARCR_07_1065,ARCR_07_1070,ARGF_01_1068,ARGF_01_1082,ARGF_01_1088,ARGF_07_1121,ARGF_07_1125,ARGF_07_1127,ARGF_07_1131,ARIG_06_1222,ARIQ_DFO_11_1071,ARPB_xx_1091,ARPB_xx_1095,ARPB_xx_1100,ARPG_19_1614,ARPI_04_1184,ARPI_08_0082,ARPI_17_1336,ARRB_16_1379,ARRB_99_1001,ARRB_99_1026,ARRB_xx_1421,ARRB_xx_1457,ARRE_02_1094,ARRE_02_1106,ARRE_06_1144,ARRE_06_1146,ARRE_06_1157,ARRE_06_1158,ARRE_06_1159,ARRE_12_1295,ARSB_09_1050,ARSB_09_1053,B96_392_SB,B96_394_SB,B96_398_SB,FMMM_RB_005_93,FMMM_RB_011_93,IMM_82_M1,MM_AB_87_005,PI_92_009_MM,RBMM01,RBMM02

done < $list
