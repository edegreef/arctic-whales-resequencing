#!/bin/bash

# Estimating ROH for each population/group separately

# adjust parameters and set# for different runs

# ROH FOR BOWHEAD
/home/degreefe/programs/plink --allow-extra-chr --vcf bowhead_RWmap_snps_rmfix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.n20.vcf.gz --set-missing-var-ids @:#\$1,\$2 --homozyg --homozyg-kb 300 --homozyg-snp 50 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 2 --homozyg-window-missing 10 --homozyg-window-threshold 0.05 --out bowhead_RWmap.n20.ROH.maf.set8 --double-id


# ROH FOR Narwhal
# divide into pops first
vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf.gz --recode --recode-INFO-all --keep EBB_group_25 --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.EBB.n25

vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf.gz --recode --recode-INFO-all --keep WBB_group_23 --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.WBB.n23

vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf.gz --recode --recode-INFO-all --keep NHB_group_9 --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.NHB.n9

/home/degreefe/programs/plink --allow-extra-chr --vcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.WBB.n23.vcf.gz --set-missing-var-ids @:#\$1,\$2 --homozyg --homozyg-kb 300 --homozyg-snp 50 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 2 --homozyg-window-missing 10 --homozyg-window-threshold 0.05 --out narwhal_WBB.ROH.maf.set8 --double-id

/home/degreefe/programs/plink --allow-extra-chr --vcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.EBB.n25.vcf.gz --set-missing-var-ids @:#\$1,\$2 --homozyg --homozyg-kb 300 --homozyg-snp 50 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 2 --homozyg-window-missing 10 --homozyg-window-threshold 0.05 --out narwhal_EBB.ROH.maf.set8 --double-id

/home/degreefe/programs/plink --allow-extra-chr --vcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.NHB.n9.vcf.gz --set-missing-var-ids @:#\$1,\$2 --homozyg --homozyg-kb 300 --homozyg-snp 50 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 2 --homozyg-window-missing 10 --homozyg-window-threshold 0.05 --out narwhal_NHB.ROH.maf.set8 --double-id
