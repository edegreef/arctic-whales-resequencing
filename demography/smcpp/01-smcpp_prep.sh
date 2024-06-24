#!/bin/bash

### preparing mask file for smcpp

# after installation finishes, then run smcpp/bin/activiate to use it
#source /home/degreefe/programs/smcpp/bin/activate


# now make vcf of filtered out sites
/home/degreefe/programs/bcftools-1.9/bcftools isec -p isec_filteredout bowhead_RWmap_allvariants.ID.vcf.gz bowhead_RWmap_snps_rmfix.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz

#isec_filteredout/0000.vcf is the one we want (snps that were filtered out)
cp isec_filteredout/0000.vcf 0000.vcf

# rename for smc input
mv 0000.vcf smcpp_removed_sites.vcf

# making bed file for masking
#export PATH=$PATH:/home/degreefe/programs/bedops/bin

# install bedops
#cd /home/degreefe/programs
#mkdir bedops
#cd bedops
#wget https://github.com/bedops/bedops/releases/download/v2.4.39/bedops_linux_x86_64-v2.4.39.tar.bz2
#tar jxvf bedops_linux_x86_64-v2.4.39.tar.bz2

export PATH=$PATH:/home/degreefe/programs/bedops/bin

# change back to your working directory

# convert to vcf to bed
switch-BEDOPS-binary-type --megarow

vcf2bed --do-not-sort --deletions < smcpp_removed_sites.vcf > bow_deletions.bed
vcf2bed --snvs < smcpp_removed_sites.vcf > bow_snps.bed

# merge bed files
bedops --everything bow_{deletions,snps}.bed > smcpp_removed_sites.bed

# sort bed
sort-bed smcpp_removed_sites.bed > smcpp_removed_sites_sorted.bed

# zip & index
bgzip smcpp_removed_sites_sorted.bed
tabix -p bed smcpp_removed_sites_sorted.bed.gz

# use this final bed file in --mask for SMC++ vcf2smc. 
