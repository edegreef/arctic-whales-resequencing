#!/bin/bash

# Script to make a vcf removing snps on X and Y chr, given list of scaffolds for each

# Input prefix of vcf and also file names for scaffolds lists for X and Y
snps=NAR_SNPS.filter1.miss.biallel.min100kb
scaf_list=final_scaffolds_Xlinked.txt
scaf_list2=final_scaffolds_Ylinked.txt

# Index vcf if not indexed already
#bgzip $snps.vcf
#tabix -p vcf $snps.vcf.gz

# Convert scaffold list to one-liner to use in bcftools
awk '{print $1}' $scaf_list | paste -s -d, - > scaf_list_line
awk '{print $1}' $scaf_list2 | paste -s -d, - > scaf_list_line2

# Set up list
list=`cat scaf_list_line`
list2=`cat scaf_list_line2`

# Filter vcf FOR these scaffolds first
bcftools filter --regions $list $snps.vcf.gz > $snps.X.vcf
bcftools filter --regions $list2 $snps.vcf.gz > $snps.Y.vcf

# Create list of CHROM and POS for the snps in the X, also removing the ##header stuff with sed
awk '{print $1, $2}' $snps.X.vcf > X_CHROM_POS
sed '/^#/d' X_CHROM_POS > X_CHROM_POS_list
rm X_CHROM_POS

awk '{print $1, $2}' $snps.Y.vcf > Y_CHROM_POS
sed '/^#/d' Y_CHROM_POS > Y_CHROM_POS_list
rm Y_CHROM_POS

# Merge these lists for one XY list to filter out
cat X_CHROM_POS_list Y_CHROM_POS_list > XY_CHROM_POS_list

# Clean up files
rm X_CHROM_POS
rm X_CHROM_POS_list
rm Y_CHROM_POS
rm Y_CHROM_POS_list

# Do the filter to remove snps on X and Y
vcftools --gzvcf $snps.vcf.gz --exclude-positions XY_CHROM_POS_list --recode --recode-INFO-all --out $snps.autosomes
mv $snps.autosomes.recode.vcf $snps.autosomes.vcf

# Zip and index
bgzip $snps.autosomes.vcf
tabix -p vcf $snps.autosomes.vcf.gz

# Can make a list of scaffolds to double check
#grep -v "^#" $snps.X.vcf | cut -f1 | sort | uniq > X_scaffold_check.txt