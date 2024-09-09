#!/bin/bash

# beagle runs into error when there are scaffolds with just 1-2 snps, need to remove these for narwhal dataset before imputing
# jump to line 20 for bowhead (did not need the extra step)

# let's make a chrom pos list of the vcf.
awk '{print $1, $2}' narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.vcf > chrom_pos
sed '/^#/d' chrom_pos > chrom_pos_nohead

# use R to pull out scaffolds with only 1-2 snps

# remove 4 scaffolds (that have 1-2 snps only)
vcftools --vcf narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.vcf --not-chr SIHG01000985.1 --not-chr SIHG01003364.1 --not-chr SIHG01003795.1 --not-chr SIHG01005513.1 --recode --recode-INFO-all --out narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm

# Zip vcf
bgzip narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.vcf
tabix -p vcf narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.vcf.gz

##### Impute! (and phase)
java -Xmx70g -jar beagle.28Jun21.220.jar gt=narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.vcf.gz iterations=20 out=narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed

# Thin snps
vcftools --gzvcf narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.vcf.gz --recode --recode-INFO-all --thin 1000 --stdout > narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000.vcf

# zip
bgzip narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000.vcf
tabix -p vcf narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000.vcf.gz

# Convert to ped/map file for LEA
/home/degreefe/programs/bcftools-1.9/bcftools view -H narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000.vcf.gz | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom-map.txt
vcftools --gzvcf narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000.vcf.gz --plink --chrom-map chrom-map.txt --out narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000