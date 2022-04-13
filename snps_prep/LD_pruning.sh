#!/bin/bash

snps=$1

# LD pruning SNPs
# --double-id b/c multiple "_" in sample name
/home/degreefe/programs/plink --vcf $snps.vcf --allow-extra-chr --indep-pairwise 50 5 0.8 --set-missing-var-ids @:#\$1,\$2 --out pruned.50window.r0.8 --double-id

# preparing the prune-in file to filter vcf
sed 's/:/\t/g' pruned.50window.r0.8.prune.in > temp.in
sed 's/...$//' temp.in > temp2.in
mv temp2.in list.pruned.50window.r0.8.prune.in
rm temp.in

# prune the snps
vcftools --vcf $snps.vcf --positions list.pruned.50window.r0.8.prune.in --recode --recode-INFO-all --out $snps.LDprunedr08

# renaming to remove the "recode" thing
mv $snps.LDprunedr08.recode.vcf $snps.LDprunedr08.vcf