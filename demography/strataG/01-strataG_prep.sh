#!/bin/bash

# example run on the command line:
# ./strataG_prep_pipeline.sh narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.n57.RB

# using vcfs that are already downsampled to population-level

# script:
snps=$1

# filter out missing data and also do the maf 0.05 filter
vcftools --gzvcf $1.vcf.gz --max-missing 1.0 --maf 0.05 --recode --recode-INFO-all --out $1.nomiss.maf05
mv $1.nomiss.maf05.recode.vcf $1.nomiss.maf05.vcf

# selecting random 25K, over 3 sets
/home/degreefe/programs/plink --thin-count 25000 --allow-extra-chr --make-bed --vcf $1.nomiss.maf05.vcf --set-missing-var-ids @:#\$1,\$2 --out $1.nomiss.maf05.thinned25K.set1 --double-id

/home/degreefe/programs/plink --thin-count 25000 --allow-extra-chr --make-bed --vcf $1.nomiss.maf05.vcf --set-missing-var-ids @:#\$1,\$2 --out $1.nomiss.maf05.thinned25K.set2 --double-id

/home/degreefe/programs/plink --thin-count 25000 --allow-extra-chr --make-bed --vcf $1.nomiss.maf05.vcf --set-missing-var-ids @:#\$1,\$2 --out $1.nomiss.maf05.thinned25K.set3 --double-id

# convert the bfiles back to vcfss
/home/degreefe/programs/plink --bfile $1.nomiss.maf05.thinned25K.set1 --recode vcf --out $1.nomiss.maf05.thinned25K.set1 --allow-extra-chr

/home/degreefe/programs/plink --bfile $1.nomiss.maf05.thinned25K.set2 --recode vcf --out $1.nomiss.maf05.thinned25K.set2 --allow-extra-chr

/home/degreefe/programs/plink --bfile $1.nomiss.maf05.thinned25K.set3 --recode vcf --out $1.nomiss.maf05.thinned25K.set3 --allow-extra-chr

# delete bfiles to clean up space
rm $1.nomiss.maf05.thinned25K.set1.bed
rm $1.nomiss.maf05.thinned25K.set1.bim
rm $1.nomiss.maf05.thinned25K.set1.fam
rm $1.nomiss.maf05.thinned25K.set1.nosex
rm $1.nomiss.maf05.thinned25K.set2.bed
rm $1.nomiss.maf05.thinned25K.set2.bim
rm $1.nomiss.maf05.thinned25K.set2.fam
rm $1.nomiss.maf05.thinned25K.set2.nosex
rm $1.nomiss.maf05.thinned25K.set3.bed
rm $1.nomiss.maf05.thinned25K.set3.bim
rm $1.nomiss.maf05.thinned25K.set3.fam
rm $1.nomiss.maf05.thinned25K.set3.nosex

