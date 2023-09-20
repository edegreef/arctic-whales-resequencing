#!/bin/bash

##################### GONE notes

#downloaded program from git: https://github.com/esrud/GONE

# Basing INPUT_PARAMETERS_FILE off Kardos et al. 2023 (https://github.com/martykardos/KillerWhaleInbreeding/tree/main/GONE) 

# GONE takes max 200 chr/scaffolds.
# Each scaffold have to be numbered 1 through max 200. 
# So if using 85 scaffolds, will have to number them 1-85.
# Can have missing data.
# Best to use plink to convert vcf to ped/map instead of vcftools b/c of the format with 0 0 0 -9


#### Narwhal
#
# starting with narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz
vcftools --gzvcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz --maf 0.05 --remove-indv 94_RAHM_IQ_162 --remove-indv ARRE_06_1164 --remove-indv B95_51_BI --remove-indv RA_HM_IQ_121 --remove-indv B96_393_SB --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.n57 
mv narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.n57.recode.vcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.n57.vcf

# Need to name each scaffold from 1-# to make format compatible with GONE.
# Narwhal has 34 scaffolds
# One way to pull out all the scaffold names:
# maybe have to remove "| sort |" between cut and uniq.
#grep -v "^#" narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.n57.vcf | cut -f1 | sort | uniq > narwhal_scaf_list_check.txt
grep -v "^#" narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.n57.vcf | cut -f1 | uniq > narwhal_scaf_list_check2.txt

# add row number to become the new scaffold names:
awk '{print $0,NR}' narwhal_scaf_list_check2.txt > narwhal_rename_scaffolds.txt

# Rename scaffolds-- in hindsight, better do to this step before separating vcfs into each pop.
/home/degreefe/programs/bcftools-1.9/bcftools annotate --rename-chrs narwhal_rename_scaffolds.txt narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.n57.vcf -Oz -o narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.n57.scafrename.vcf.gz

# Convert to ped/map format (*important to do this with plink and not vcftools)
/home/degreefe/programs/plink --vcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.n57.scafrename.vcf.gz --recode --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.n57.scafrename --double-id --chr-set 34

# narwhal with maf filter (n57): 2,207,956 snps
# narwhal without maf filter (n57): 4,962,230 snps


############# Bowhead has more than 200 scaffolds so need extra step.

#### Bowhead whale
# starting with bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz. doing maf 0.05 and removing the two kin pairs
vcftools --gzvcf bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz --maf 0.05 --remove-indv ARBMGH_2002_001 --remove-indv 99_01 --recode --recode-INFO-all --out bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.n21 
mv bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.n21.recode.vcf bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.n21.vcf

# Need to name each scaffold from 1-# to make format compatible with GONE.
# Before renaming scaffolds, select top 200 scaffolds
# Run the bowhead_extract_top200scaf.sh script
# Use difference scaffold list here if want to use different 200 scaffolds.
./bowhead_extract_200scaf.sh

## for selecting a random 200 scaffold set, start with the bowhead_min100kb.txt list, then run
#shuf -n 200 bowhead_min100kb.txt > bowhead_random200scaf_set2.txt
#shuf -n 200 bowhead_min100kb.txt > bowhead_random200scaf_set3.txt
## and then edit the ./bowhead_extract_200scaf.sh script


# Do a quick check on the scaffolds
grep -v "^#" bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.n21.top200scaf.vcf | cut -f1 | uniq > bowhead_scaf_list_check.txt

# add row number to become the new scaffold names:
awk '{print $0,NR}' bowhead_scaf_list_check.txt > bowhead_rename_scaffolds.txt

# For situation with many scaffolds, temporarly add a "s" for each number (the --allow-extra-chr in plink doesn't work for raw numbers...) will remove this s later because then GONE only takes number values ... ugh.
# Would be fine if plink doesn't max number chr at 95 -- specific to bowhead maybe?
awk '$2="s"$2' bowhead_rename_scaffolds.txt > bowhead_rename_scaffolds_s.txt


# Rename scaffolds
# top 200
/home/degreefe/programs/bcftools-1.9/bcftools annotate --rename-chrs bowhead_rename_scaffolds_s.txt bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.n21.top200.vcf -Oz -o bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.n21.top200.scafrename.vcf.gz

# Convert to ped/map format (*important to do this with plink and not vcftools)
# using allow-extra-chr instead of --chr-set (b/c of the 95 cap)
/home/degreefe/programs/plink --vcf bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.n21.top200.scafrename.vcf.gz --recode --out bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.n21.top200.scafrename --double-id --allow-extra-chr

# now need to remove the "s" that we added to bypass the 95 cap.
# removing the first letter of every line in map file
cut -c2- bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.n21.top200.scafrename.map > bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.n21.top200.scafrename.temp.map

# yay it worked
# clean up a bit
rm bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.n21.top200.scafrename.map
mv bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.n21.top200.scafrename.temp.map  bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.n21.top200scaf.map

