#!/bin/bash

# Looking at snp statistics

# Vcf prefix names for input and output
all_variants=BOW_allvariantcalls.ID
snps=BOW_SNPS

# Filter out indels 
vcftools --gzvcf $all_variants.vcf.gz --remove-indels --recode --recode-INFO-all --out $snps
mv $snps.recode.vcf $snps.vcf

# Zip and index vcf if not done already
bgzip $snps.vcf
tabix -p vcf $snps.vcf.gz

# Get snp info using bcftools 
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/MQ\t%INFO/QD\t%INFO/TC\t%INFO/FR\n' $snps.vcf.gz -o $snps.info

# Allele frequency (outputs .frq file)
vcftools --gzvcf $snps.vcf.gz --out $snps --freq2 --max-alleles 2

# Proportion of missing data per individual (outputs .imiss file)
vcftools --gzvcf $snps.vcf.gz --out $snps --missing-indv 

# Propotion of missing data per site (outputs .lmiss file)
vcftools --gzvcf $snps.vcf.gz --out $snps --missing-site

# Assess sites for Hardy-Weinberg Equilibrium (outputs .hwe file)
vcftools --gzvcf $snps.vcf.gz --out $snps --hardy