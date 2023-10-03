#!/bin/bash

# Plotting roh using outputs from PLINK
# Ran ROH for each population/group separately


# Filter vcf for each pop/region (probably a more efficient way to do this, but here goes)
vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.vcf.gz --keep indiv_AB --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.AB

vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.vcf.gz --keep indiv_BI --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.BI

vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.vcf.gz --keep indiv_CR --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.CR

vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.vcf.gz --keep indiv_GF --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.GF

vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.vcf.gz --keep indiv_IG --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.IG

vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.vcf.gz --keep indiv_PB --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.PB

vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.vcf.gz --keep indiv_PG --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.PG

vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.vcf.gz --keep indiv_PI --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.PI

vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.vcf.gz --keep indiv_RB --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.RB

vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.vcf.gz --keep indiv_RE --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.RE

vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.vcf.gz --keep indiv_SB --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.SB

vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.vcf.gz --keep indiv_SB --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.SB

vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.vcf.gz --keep indiv_GF --keep indiv_RE --keep indiv_SB --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.NOR

vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.vcf.gz --keep indiv_PB --keep indiv_IG --keep indiv_AB --keep indiv_PG --keep indiv_PI --keep indiv_CR --keep indiv_BI --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.MID

# ROH line for just one site/group:
#/home/degreefe/programs/plink --allow-extra-chr --vcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.AB.vcf --set-missing-var-ids @:#\$1,\$2 --homozyg --homozyg-kb 50 --homozyg-snp 50 --out AB --double-id

# Run ROH, TAJD, and PI for all the sites/groups using list of site abbrevs in "site_list"
while read site
do

# ROH (runs of homozygosity)
/home/degreefe/programs/plink --allow-extra-chr --vcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.$site.vcf.gz --set-missing-var-ids @:#\$1,\$2 --homozyg --homozyg-kb 300 --homozyg-snp 50 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 3 --homozyg-window-missing 10 --homozyg-window-threshold 0.05 --out $site --double-id

# TAJD (Tajima's D)
vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.$site.vcf.gz --TajimaD 10000 --out $site.10kb

# PI (Nucleotide diversity)
vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.$site.vcf.gz --window-pi 10000 --out $site.10kb

done < site_list