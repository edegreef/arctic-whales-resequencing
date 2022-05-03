#!/bin/bash

### preparing mask file for smcpp

# Notes on installing SMC++ and preparing input files

# download smc++
wget https://github.com/popgenmethods/smcpp/releases/download/v1.15.2/smcpp-1.15.2-Linux-x86_64.sh

chmod +x smcpp-1.15.2-Linux-x86_64.sh

# install
./smcpp-1.15.2-Linux-x86_64.sh

# go through reading and acceptng license - answer 'yes' at the end, then program will be installed under /home/edegreef/smcpp

# after installation finishes, then run smcpp/bin/activiate to use it
source /home/edegreef/smcpp/bin/activate

# check it works
smc++ {options}


# prep input files
# zip and index if not already done (for both RAW variants, and filtered variants)
#bgzip BOW_SNPS.filter1.miss.biallel.min100kb.autosomes.vcf
#tabix -p vcf BOW_SNPS.filter1.miss.biallel.min100kb.autosomes.vcf.gz

cd /home/degreefe/whales/BOW_snps/smcpp

# now make vcf of filtered out sites
/home/degreefe/programs/bcftools-1.9/bcftools isec -p isec_filteredout BOW_allvariantcalls.ID.vcf.gz BOW_SNPS.filter1.miss.biallel.min100kb.autosomes.vcf.gz

# isec_filteredout/0000.vcf is the one we want (snps that were filtered out)
cp isec_filteredout/0000.vcf 0000.vcf

# rename for smc input
mv 0000.vcf smcpp_removed_sites.vcf

# making bed file for masking
export PATH=$PATH:/home/degreefe/programs/bedops/bin

# install bedops
cd /home/degreefe/programs
mkdir bedops
cd bedops
wget https://github.com/bedops/bedops/releases/download/v2.4.39/bedops_linux_x86_64-v2.4.39.tar.bz2
tar jxvf bedops_linux_x86_64-v2.4.39.tar.bz2

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
