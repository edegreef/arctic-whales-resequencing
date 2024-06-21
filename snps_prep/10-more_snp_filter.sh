#!/bin/bash

# Converting to bfiles, running kinship, and filtering out close kin/dups/high miss, then convert to bfile again, also convert to map/ped for LEA. A few lines specific to bowhead/narwhal and need to be adjusted depending which one is running.

snps=bowhead_RWmap_snps_rmfix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08
#narwhal: narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08
#bowhead: bowhead_RWmap_snps_rmfix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08

n=n20
#narwhal: n57 
#bowhead: n20


# Kinship
# Convert to bim/bed/fam
/home/degreefe/programs/plink --allow-extra-chr --make-bed --vcf $snps.vcf.gz --set-missing-var-ids @:#\$1,\$2 --out $snps --double-id

# Run plink pairwise estimates, will use the .genome and .imiss files for plinkQC in R
/home/degreefe/programs/plink --allow-extra-chr --bfile $snps  --genome --out $snps 
/home/degreefe/programs/plink --allow-extra-chr --bfile $snps  --missing --out $snps 

# BOWHEAD: remove sample 99_01 b/c close kin
vcftools --gzvcf $snps.vcf.gz --remove-indv 99_01 --recode --recode-INFO-all --out $snps.$n

# NARWHAL: remove dups, close kin, and high miss
#vcftools --gzvcf $snps.vcf.gz --remove-indv 94_RAHM_IQ_162 --remove-indv ARRE_06_1164 --remove-indv B95_51_BI --remove-indv B96_393_SB --remove-indv RA_HM_IQ_121 --recode --recode-INFO-all --out $snps.$n

# Convert to bim/bed/fam
/home/degreefe/programs/plink --allow-extra-chr --make-bed --vcf $snps.$n.vcf  --set-missing-var-ids @:#\$1,\$2 --out $snps.$n --double-id

# Prep for LEA - Convert to ped/map 
/home/degreefe/programs/bcftools-1.9/bcftools view -H $snps.$n.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom-map.txt
vcftools --vcf $snps.$n.vcf --plink --chrom-map chrom-map.txt --out $snps.$n