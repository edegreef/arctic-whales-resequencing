#!/bin/bash

####### 1) Narwhal
# Remove kin/duplicates
vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.vcf.gz --remove-indv 94_RAHM_IQ_162 --remove-indv RA_HM_IQ_121 --remove-indv B96_393_SB --remove-indv ARRE_06_1164 --remove-indv B95_51_BI --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.n57

# Minor allele count filter
vcftools --vcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.n57.vcf --mac 2 --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.n57.mac2

# 3,777,422 out of 5,002,361 snps

# Remove SNPs with more than 10% missing data (max missingness 0.1, but here in vcftools is flipped (1=no missing, so 0.9 for max-missing of 10%)
vcftools --vcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.n57.mac2.vcf --max-missing 0.9 --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.n57.mac2.miss01

# 3,761,357 snps

# Zip and index vcf
bgzip narwhal_snps.filter1.miss.biallel.min100kb.autosomes.n57.mac2.miss01.vcf
tabix -p vcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.n57.mac2.miss01.vcf.gz

# Impute (and phase) with beagle
java -Xmx70g -jar beagle.28Jun21.220.jar gt=narwhal_snps.filter1.miss.biallel.min100kb.autosomes.n57.mac2.miss01.vcf.gz iterations=20 out=narwhal_snps.filter1.miss.biallel.min100kb.autosomes.n57.mac2.miss01.imputed




####### 2) Bowhead whale
# Remove two kin pairs
vcftools --gzvcf bowhead_snps.filter1.miss.biallel.min100kb.autosomes.vcf.gz --remove-indv 99_01 --remove-indv ARBMGH_2002_001 --recode --recode-INFO-all --out bowhead_snps.filter1.miss.biallel.min100kb.autosomes.n21

# Minor allele count filter
vcftools --vcf bowhead_snps.filter1.miss.biallel.min100kb.autosomes.n21.vcf --mac 2 --recode --recode-INFO-all --out bowhead_snps.filter1.miss.biallel.min100kb.autosomes.n21.mac2

# 8,078,637 out of 9,848,912 snps

# Max missingness 10%
vcftools --vcf bowhead_snps.filter1.miss.biallel.min100kb.autosomes.n21.mac2.vcf --max-missing 0.9 --recode --recode-INFO-all --out bowhead_snps.filter1.miss.biallel.min100kb.autosomes.n21.mac2.miss01

# 7,999,069 snps

# Zip and index vcf
bgzip bowhead_snps.filter1.miss.biallel.min100kb.autosomes.n21.mac2.miss01.vcf
tabix -p vcf bowhead_snps.filter1.miss.biallel.min100kb.autosomes.n21.mac2.miss01.vcf.gz

# Impute (and phase) with beagle
java -Xmx70g -jar beagle.28Jun21.220.jar gt=bowhead_snps.filter1.miss.biallel.min100kb.autosomes.n21.mac2.miss01.vcf.gz iterations=20 out=bowhead_snps.filter1.miss.biallel.min100kb.autosomes.n21.mac2.miss01.imputed